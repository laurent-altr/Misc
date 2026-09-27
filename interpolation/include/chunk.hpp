// Chunked interpolation: m independent lanes, each with its own x and its own
// current interval k (the initial guess, updated in place). Two loops, both
// written so that the compiler can vectorize them:
//   update(): moves each k[i] by one interval per pass until every lane is in
//             place (almost always 0 or 1 pass when the guess is good);
//   eval_*(): gathers the data of interval k[i] and interpolates.
#pragma once
#include <vector>

#include "interp.hpp"

namespace chunk {

// Table in structure-of-arrays layout: knots, plus per-interval data for the
// compensated formula (same values as interp::make_comp_seg).
template <class T> struct Table {
    std::vector<T> x, y;             // knots
    std::vector<T> s, s_lo, sh, sl;  // per interval: slope, its correction, Veltkamp split of the slope
    std::vector<T> slope;            // per interval: plain slope, for the naive formula with a precomputed slope
    int intervals() const { return int(x.size()) - 1; }
};

template <class T> Table<T> make_table(const std::vector<T>& X, const std::vector<T>& Y) {
    Table<T> t{X, Y, {}, {}, {}, {}, {}};
    for (std::size_t k = 0; k + 1 < X.size(); ++k) {
        const auto g = interp::make_comp_seg(X[k], X[k + 1], Y[k], Y[k + 1]);
        t.s.push_back(g.s);
        t.s_lo.push_back(g.s_lo);
        t.sh.push_back(g.s_hi_half);
        t.sl.push_back(g.s_lo_half);
        t.slope.push_back((Y[k + 1] - Y[k]) / (X[k + 1] - X[k]));
    }
    return t;
}

// Moves each k[i] until X[k] <= x[i] < X[k+1] (k stays within [0, intervals-1]).
template <class T> void update(const Table<T>& t, const T* __restrict x, int* __restrict k, int m) {
    const T* __restrict X = t.x.data();
    const int last = t.intervals() - 1;
    for (;;) {
        int moved = 0;
        for (int i = 0; i < m; ++i) {
            const int ki = k[i];
            const int up = (ki < last) & (x[i] >= X[ki + 1]);
            const int down = (ki > 0) & (x[i] < X[ki]);
            k[i] = ki + up - down;
            moved |= up | down;
        }
        if (!moved) return;
    }
}

// Naive formula anchored on the left knot: y0 + (x - x0) * dy / dx.
template <class T> void eval_naive(const Table<T>& t, const T* __restrict x, const int* __restrict k, T* __restrict y, int m) {
    const T* __restrict X = t.x.data();
    const T* __restrict Y = t.y.data();
    for (int i = 0; i < m; ++i) {
        const int j = k[i];
        y[i] = Y[j] + (x[i] - X[j]) * (Y[j + 1] - Y[j]) / (X[j + 1] - X[j]);
    }
}

// Naive formula with the slope precomputed: y0 + (x - x0) * slope.
template <class T> void eval_naive_slope(const Table<T>& t, const T* __restrict x, const int* __restrict k, T* __restrict y, int m) {
    const T* __restrict X = t.x.data();
    const T* __restrict Y = t.y.data();
    const T* __restrict S = t.slope.data();
    for (int i = 0; i < m; ++i) {
        const int j = k[i];
        y[i] = Y[j] + (x[i] - X[j]) * S[j];
    }
}

// Compensated nearest-bound formula (H' with precomputed interval data).
// Same bits as interp::eval_comp_seg.
// The loop bodies are written out in full on purpose: with gcc 13, calling
// an inlined per-lane helper here prevents vectorization ("no vectype").
template <class T>
INTERP_NO_CONTRACT void eval_comp(const Table<T>& t, const T* __restrict x, const int* __restrict k, T* __restrict y, int m) {
    const T* __restrict X = t.x.data();
    const T* __restrict Y = t.y.data();
    const T* __restrict S = t.s.data();
    const T* __restrict SL = t.s_lo.data();
    const T* __restrict SH = t.sh.data();
    const T* __restrict SS = t.sl.data();
    for (int i = 0; i < m; ++i) {
        const int j = k[i];
        const T xi = x[i], x0 = X[j], x1 = X[j + 1];
        // Nearest knot as an index (j or j + 1): one gather, no conditional
        // load (a conditional load also prevents vectorization).
        const int ja = j + int(!(xi - x0 < x1 - xi));
        const T xa = X[ja], ya = Y[ja];
        T h, h_e, hh, hl, sum, sum_e;
        interp::two_sum(xi, -xa, h, h_e);
        const T p = h * S[j];
        interp::split(h, hh, hl);
        const T p_e = ((hh * SH[j] - p) + hh * SS[j] + hl * SH[j]) + hl * SS[j];
        interp::two_sum(ya, p, sum, sum_e);
        y[i] = sum + (sum_e + (p_e + (h * SL[j] + h_e * S[j])));
    }
}

// Compensated formula anchored on the left knot (L): same accuracy as
// eval_comp (the anchor does not matter once every rounding error is
// compensated), without the nearest-knot choice: no compare, two fewer gathers.
// Same bits as interp::left_comp_nofma.
template <class T>
INTERP_NO_CONTRACT void eval_comp_left(const Table<T>& t, const T* __restrict x, const int* __restrict k, T* __restrict y, int m) {
    const T* __restrict X = t.x.data();
    const T* __restrict Y = t.y.data();
    const T* __restrict S = t.s.data();
    const T* __restrict SL = t.s_lo.data();
    const T* __restrict SH = t.sh.data();
    const T* __restrict SS = t.sl.data();
    for (int i = 0; i < m; ++i) {
        const int j = k[i];
        T h, h_e, hh, hl, sum, sum_e;
        interp::two_sum(x[i], -X[j], h, h_e);
        const T p = h * S[j];
        interp::split(h, hh, hl);
        const T p_e = ((hh * SH[j] - p) + hh * SS[j] + hl * SH[j]) + hl * SS[j];
        interp::two_sum(Y[j], p, sum, sum_e);
        y[i] = sum + (sum_e + (p_e + (h * SL[j] + h_e * S[j])));
    }
}

}  // namespace chunk
