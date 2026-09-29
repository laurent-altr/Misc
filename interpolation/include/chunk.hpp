// Chunked interpolation: m independent lanes, each with its own x and its own
// current interval k (the initial guess, updated in place). Two loops, both
// written so that the compiler can vectorize them:
//   update(): moves each k[i] by one interval per pass until every lane is in
//             place (almost always 0 or 1 pass when the guess is good);
//   eval_*(): gathers the data of interval k[i] and interpolates.
#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
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
    INTERP_STRICT_FP
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
    INTERP_STRICT_FP
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

#if defined(__FMA__)
// L with hardware fma for the exact product (2 instructions instead of the
// Veltkamp/Dekker product). Same bits as eval_comp_left.
template <class T>
INTERP_NO_CONTRACT void eval_comp_left_fma(const Table<T>& t, const T* __restrict x, const int* __restrict k, T* __restrict y, int m) {
    INTERP_STRICT_FP
    const T* __restrict X = t.x.data();
    const T* __restrict Y = t.y.data();
    const T* __restrict S = t.s.data();
    const T* __restrict SL = t.s_lo.data();
    for (int i = 0; i < m; ++i) {
        const int j = k[i];
        T h, h_e, sum, sum_e;
        interp::two_sum(x[i], -X[j], h, h_e);
        const T p = h * S[j];
        const T p_e = std::fma(h, S[j], -p);
        interp::two_sum(Y[j], p, sum, sum_e);
        y[i] = sum + (sum_e + (p_e + (h * SL[j] + h_e * S[j])));
    }
}
#endif

// Per-lane cache of the current interval, in structure-of-arrays layout, for
// x that changes slowly: evaluation then reads contiguous memory (no gather),
// and the interval check is a contiguous compare with an early exit. Only
// lanes whose x left their interval are refreshed (scalar).
template <class T, int M> struct Lanes {
    alignas(64) T x0[M];      // anchor (left knot)
    alignas(64) T lo[M];      // interval bounds for the check (lowest / max at the table ends,
    alignas(64) T hi[M];      //   so lanes outside the table extrapolate the end intervals)
    alignas(64) T y0[M];
    alignas(64) T s[M];       // slope
    alignas(64) T s_lo[M];    // slope correction
    alignas(64) T sh[M];      // Veltkamp split of the slope (fma-free product)
    alignas(64) T sl[M];
    alignas(64) int k[M];     // current interval
};

template <class T, int M> void fill_lane(const Table<T>& t, Lanes<T, M>& c, int i, int j) {
    // Finite sentinels rather than infinities, which are undefined under
    // -ffast-math / icpx -fp-model=fast.
    c.k[i] = j;
    c.x0[i] = t.x[j];
    c.lo[i] = j == 0 ? std::numeric_limits<T>::lowest() : t.x[j];
    c.hi[i] = j == t.intervals() - 1 ? std::numeric_limits<T>::max() : t.x[j + 1];
    c.y0[i] = t.y[j];
    c.s[i] = t.s[j];
    c.s_lo[i] = t.s_lo[j];
    c.sh[i] = t.sh[j];
    c.sl[i] = t.sl[j];
}

template <class T> int hunt(const Table<T>& t, T x, int k) {
    const T* X = t.x.data();
    const int last = t.intervals() - 1;
    while (k > 0 && x < X[k]) --k;
    while (k < last && x >= X[k + 1]) ++k;
    return k;
}

// Initial fill: k[i] is the initial guess.
template <class T, int M> void init(const Table<T>& t, Lanes<T, M>& c, const T* x, const int* k) {
    for (int i = 0; i < M; ++i) fill_lane(t, c, i, hunt(t, x[i], k[i]));
}

// Returns the number of lanes that changed interval.
template <class T, int M> int refresh(const Table<T>& t, Lanes<T, M>& c, const T* __restrict x) {
    // The flag uses an integer as wide as T: mixing 32-bit ints with double
    // comparisons prevents vectorization with gcc 13.
    using I = std::conditional_t<sizeof(T) == 4, std::int32_t, std::int64_t>;
    I any = 0;
    for (int i = 0; i < M; ++i) any |= I(x[i] < c.lo[i]) | I(x[i] >= c.hi[i]);
    if (!any) return 0;  // early exit: the common case when x changes slowly
    int moved = 0;
    for (int i = 0; i < M; ++i) {
        if (x[i] < c.lo[i] || x[i] >= c.hi[i]) {
            fill_lane(t, c, i, hunt(t, x[i], c.k[i]));
            ++moved;
        }
    }
    return moved;
}

// L from the cache, fma-free. Same bits as eval_comp_left.
template <class T, int M>
INTERP_NO_CONTRACT void eval_cached(const Lanes<T, M>& c, const T* __restrict x, T* __restrict y) {
    INTERP_STRICT_FP
    for (int i = 0; i < M; ++i) {
        T h, h_e, hh, hl, sum, sum_e;
        interp::two_sum(x[i], -c.x0[i], h, h_e);
        const T p = h * c.s[i];
        interp::split(h, hh, hl);
        const T p_e = ((hh * c.sh[i] - p) + hh * c.sl[i] + hl * c.sh[i]) + hl * c.sl[i];
        interp::two_sum(c.y0[i], p, sum, sum_e);
        y[i] = sum + (sum_e + (p_e + (h * c.s_lo[i] + h_e * c.s[i])));
    }
}

#if defined(__FMA__)
// L from the cache, with hardware fma. Same bits.
template <class T, int M>
INTERP_NO_CONTRACT void eval_cached_fma(const Lanes<T, M>& c, const T* __restrict x, T* __restrict y) {
    INTERP_STRICT_FP
    for (int i = 0; i < M; ++i) {
        T h, h_e, sum, sum_e;
        interp::two_sum(x[i], -c.x0[i], h, h_e);
        const T p = h * c.s[i];
        const T p_e = std::fma(h, c.s[i], -p);
        interp::two_sum(c.y0[i], p, sum, sum_e);
        y[i] = sum + (sum_e + (p_e + (h * c.s_lo[i] + h_e * c.s[i])));
    }
}
#endif

// Naive formula from the cache (for a fair comparison): y0 + (x - x0) * slope.
// Uses the rounded slope s, which equals the plain precomputed slope.
template <class T, int M> void eval_cached_naive(const Lanes<T, M>& c, const T* __restrict x, T* __restrict y) {
    for (int i = 0; i < M; ++i) y[i] = c.y0[i] + (x[i] - c.x0[i]) * c.s[i];
}

}  // namespace chunk
