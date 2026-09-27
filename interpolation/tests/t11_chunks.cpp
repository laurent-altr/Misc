// Chunks of 128 independent lanes, each with its own x and its own current
// interval (a good initial guess, updated in place): the pattern of the target
// application. Compares, in ns per point:
//   scalar : one lane at a time, walk from the previous interval then interpolate
//   chunked: update all 128 intervals, then interpolate all 128 lanes
//            (two loops the compiler can vectorize)
// for the naive formula, the naive formula with a precomputed slope, and the
// compensated formula H' with precomputed interval data.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>

#include "chunk.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "report.hpp"

constexpr int M = 128;   // lanes per chunk
constexpr int C = 512;   // chunks per run

template <class Fn> double best_ns(Fn&& fn, long points, int reps = 60) {
    double best = 1e30;
    for (int r = 0; r < reps; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        fn();
        auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count() / points);
    }
    return best;
}

template <class T> std::size_t hunt(const T* X, int nint, T x, std::size_t k) {
    while (k > 0 && x < X[k]) --k;
    while (k + 1 < std::size_t(nint) && x >= X[k + 1]) ++k;
    return k;
}

// One lane at a time (scalar), as in t10's hunt pattern.
template <class T>
__attribute__((noinline)) void scalar_naive(const chunk::Table<T>& t, const T* __restrict x, int* __restrict k, T* __restrict y) {
    const T* X = t.x.data();
    const T* Y = t.y.data();
    for (int i = 0; i < M; ++i) {
        const std::size_t j = k[i] = int(hunt(X, t.intervals(), x[i], k[i]));
        y[i] = interp::left(X[j], X[j + 1], Y[j], Y[j + 1], x[i]);
    }
}
template <class T>
__attribute__((noinline)) void scalar_comp(const chunk::Table<T>& t, const std::vector<interp::CompSeg<T>>& segs,
                                           const T* __restrict x, int* __restrict k, T* __restrict y) {
    for (int i = 0; i < M; ++i) {
        const std::size_t j = k[i] = int(hunt(t.x.data(), t.intervals(), x[i], k[i]));
        y[i] = interp::eval_comp_seg(segs[j], x[i]);
    }
}

// Chunked: update, then evaluate.
template <class T> __attribute__((noinline)) void do_update(const chunk::Table<T>& t, const T* x, int* k) { chunk::update(t, x, k, M); }
template <class T> __attribute__((noinline)) void do_naive(const chunk::Table<T>& t, const T* x, const int* k, T* y) { chunk::eval_naive(t, x, k, y, M); }
template <class T> __attribute__((noinline)) void do_slope(const chunk::Table<T>& t, const T* x, const int* k, T* y) { chunk::eval_naive_slope(t, x, k, y, M); }
template <class T> __attribute__((noinline)) void do_comp(const chunk::Table<T>& t, const T* x, const int* k, T* y) { chunk::eval_comp(t, x, k, y, M); }

template <class T> void run(report::Checks& checks) {
    std::printf("\n=== %s ===\n", fp::type_name<T>());
    std::mt19937_64 g(21);
    std::uniform_real_distribution<double> u(0, 1), v(-100, 100);
    std::vector<T> X, Y;
    { double x = 0; for (int k = 0; k < 1024; ++k) { X.push_back(T(x)); Y.push_back(T(v(g))); x += 0.01 + u(g); } }
    const auto tab = chunk::make_table(X, Y);
    std::vector<interp::CompSeg<T>> segs;
    for (int k = 0; k < tab.intervals(); ++k) segs.push_back(interp::make_comp_seg(X[k], X[k + 1], Y[k], Y[k + 1]));

    // Each lane walks independently: steps of up to 5 % of the mean knot spacing.
    const double span = double(X.back()), mean = span / (X.size() - 1);
    std::vector<T> xs(std::size_t(C) * M);
    std::vector<int> k0(M);
    for (int i = 0; i < M; ++i) {
        double x = u(g) * span * 0.999;
        k0[i] = int(std::upper_bound(X.begin(), X.end(), T(x)) - X.begin()) - 1;
        for (int c = 0; c < C; ++c) {
            x = std::clamp(x + (u(g) - 0.5) * 0.1 * mean, 0.0, span * 0.999);
            xs[std::size_t(c) * M + i] = T(x);
        }
    }
    std::vector<T> y(std::size_t(C) * M), y2(std::size_t(C) * M);
    std::vector<int> k(M);
    const long P = long(C) * M;

    // Correctness: the chunked compensated result equals eval_comp_seg bit for bit,
    // and every lane ends in the interval containing its x.
    {
        k = k0;
        long bad_k = 0, bad_y = 0, moves = 0, passes_gt1 = 0;
        for (int c = 0; c < C; ++c) {
            const T* xc = &xs[std::size_t(c) * M];
            std::vector<int> before = k;
            chunk::update(tab, xc, k.data(), M);
            chunk::eval_comp(tab, xc, k.data(), &y[std::size_t(c) * M], M);
            for (int i = 0; i < M; ++i) {
                bad_k += !(X[k[i]] <= xc[i] && xc[i] < X[k[i] + 1]);
                const T ref = interp::eval_comp_seg(segs[k[i]], xc[i]);
                bad_y += fp::ordered(ref) != fp::ordered(y[std::size_t(c) * M + i]);
                moves += k[i] != before[i];
                passes_gt1 += std::abs(k[i] - before[i]) > 1;
            }
        }
        std::printf("lanes changing interval per chunk: %.1f %% (%.3f %% by more than one)\n", 100.0 * moves / P,
                    100.0 * passes_gt1 / P);
        checks.expect(bad_k == 0, std::string("update() finds the right interval (") + fp::type_name<T>() + ")");
        checks.expect(bad_y == 0, std::string("chunked H' identical to eval_comp_seg (") + fp::type_name<T>() + ")");

    }

    auto time = [&](auto body) {
        return best_ns([&] {
            k = k0;
            for (int c = 0; c < C; ++c) body(&xs[std::size_t(c) * M], k.data(), &y[std::size_t(c) * M]);
        }, P);
    };
    const double upd = time([&](const T* x, int* kk, T*) { do_update(tab, x, kk); });
    const double s_naive = time([&](const T* x, int* kk, T* yy) { scalar_naive(tab, x, kk, yy); });
    const double s_comp = time([&](const T* x, int* kk, T* yy) { scalar_comp(tab, segs, x, kk, yy); });
    const double c_naive = time([&](const T* x, int* kk, T* yy) { do_update(tab, x, kk); do_naive(tab, x, kk, yy); });
    const double c_slope = time([&](const T* x, int* kk, T* yy) { do_update(tab, x, kk); do_slope(tab, x, kk, yy); });
    const double c_comp = time([&](const T* x, int* kk, T* yy) { do_update(tab, x, kk); do_comp(tab, x, kk, yy); });

    report::Table t({"formula", "scalar", "chunked", "(update)"}, 12);
    report::Csv csv(std::string("t11_chunks_") + fp::type_name<T>(), "formula,scalar_ns,chunked_ns,update_ns");
    auto f = [](double v) { return v > 0 ? report::num(v, "%.2f") : std::string("-"); };
    t.line({"naive", f(s_naive), f(c_naive), f(upd)});
    t.line({"naive, slope", "-", f(c_slope), f(upd)});
    t.line({"H' (pre)", f(s_comp), f(c_comp), f(upd)});
    csv.row("naive", s_naive, c_naive, upd);
    csv.row("naive_slope", 0, c_slope, upd);
    csv.row("comp_pre", s_comp, c_comp, upd);
    std::printf("(ns per point, %d chunks of %d lanes)\n", C, M);
}

int main() {
    std::printf("t11: chunks of %d lanes with per-lane interval guesses\n", M);
    report::Checks checks;
    run<float>(checks);
    run<double>(checks);
    return checks.exit_code();
}
