// Optimizing L for chunks of 128 lanes with slowly changing x:
//   table : update(k) + evaluation gathering the interval data from the table
//   cached: per-lane cache of the interval data (contiguous loads), refreshed
//           only for the lanes whose x left their interval (early exit)
// each fma-free and with hardware fma (when compiled with FMA support), and
// the naive formula with a precomputed slope under the same two layouts.
// All L variants must give the same bits as interp::left_comp_nofma.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>

#include "chunk.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "report.hpp"

constexpr int M = 128;
constexpr int C = 512;

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

template <class T> using L = chunk::Lanes<T, M>;
#define NOINLINE __attribute__((noinline))
template <class T> NOINLINE void t_update(const chunk::Table<T>& t, const T* x, int* k) { chunk::update(t, x, k, M); }
template <class T> NOINLINE void t_left(const chunk::Table<T>& t, const T* x, const int* k, T* y) { chunk::eval_comp_left(t, x, k, y, M); }
template <class T> NOINLINE void t_slope(const chunk::Table<T>& t, const T* x, const int* k, T* y) { chunk::eval_naive_slope(t, x, k, y, M); }
template <class T> NOINLINE int c_refresh(const chunk::Table<T>& t, L<T>& c, const T* x) { return chunk::refresh(t, c, x); }
template <class T> NOINLINE void c_left(const L<T>& c, const T* x, T* y) { chunk::eval_cached(c, x, y); }
template <class T> NOINLINE void c_naive(const L<T>& c, const T* x, T* y) { chunk::eval_cached_naive(c, x, y); }
#if defined(__FMA__)
template <class T> NOINLINE void t_left_fma(const chunk::Table<T>& t, const T* x, const int* k, T* y) { chunk::eval_comp_left_fma(t, x, k, y, M); }
template <class T> NOINLINE void c_left_fma(const L<T>& c, const T* x, T* y) { chunk::eval_cached_fma(c, x, y); }
#endif

template <class T> void run(report::Checks& checks, double step, const char* label) {
    std::printf("\n=== %s, %s x (steps up to %.1f %% of the knot spacing) ===\n", fp::type_name<T>(), label, 100 * step);
    std::mt19937_64 g(21);
    std::uniform_real_distribution<double> u(0, 1), v(-100, 100);
    std::vector<T> X, Y;
    { double x = 0; for (int k = 0; k < 1024; ++k) { X.push_back(T(x)); Y.push_back(T(v(g))); x += 0.01 + u(g); } }
    const auto tab = chunk::make_table(X, Y);
    const double span = double(X.back()), mean = span / (X.size() - 1);
    std::vector<T> xs(std::size_t(C) * M);
    std::vector<int> k0(M);
    for (int i = 0; i < M; ++i) {
        double x = u(g) * span * 0.999;
        k0[i] = int(std::upper_bound(X.begin(), X.end(), T(x)) - X.begin()) - 1;
        for (int c = 0; c < C; ++c) {
            x = std::clamp(x + (u(g) - 0.5) * 2 * step * mean, 0.0, span * 0.999);
            xs[std::size_t(c) * M + i] = T(x);
        }
    }
    const long P = long(C) * M;
    std::vector<T> y(P), ref(P);
    std::vector<int> k(M);
    auto* cache = new L<T>;
    auto chunk_x = [&](int c) { return &xs[std::size_t(c) * M]; };
    auto chunk_y = [&](std::vector<T>& v, int c) { return &v[std::size_t(c) * M]; };

    // Reference: left_comp_nofma on the right interval.
    for (int c = 0; c < C; ++c)
        for (int i = 0; i < M; ++i) {
            const T xi = chunk_x(c)[i];
            const int j = std::clamp(int(std::upper_bound(X.begin(), X.end(), xi) - X.begin()) - 1, 0, tab.intervals() - 1);
            ref[std::size_t(c) * M + i] = interp::left_comp_nofma(X[j], X[j + 1], Y[j], Y[j + 1], xi);
        }
    auto same = [&](const char* what) {
        long bad = 0;
        for (long i = 0; i < P; ++i) bad += fp::ordered(y[i]) != fp::ordered(ref[i]);
        checks.expect(bad == 0, std::string(what) + " identical to left_comp_nofma (" + fp::type_name<T>() + ", " + label + ")");
    };

    // Runners (also used for the correctness checks).
    auto run_table = [&](auto eval) {
        k = k0;
        for (int c = 0; c < C; ++c) { t_update(tab, chunk_x(c), k.data()); eval(chunk_x(c), k.data(), chunk_y(y, c)); }
    };
    long moved = 0;
    auto run_cached = [&](auto eval) {
        chunk::init(tab, *cache, chunk_x(0), k0.data());
        moved = 0;
        for (int c = 0; c < C; ++c) { moved += c_refresh(tab, *cache, chunk_x(c)); eval(chunk_x(c), chunk_y(y, c)); }
    };

    auto tl = [&](const T* x, const int* kk, T* yy) { t_left(tab, x, kk, yy); };
    auto ts = [&](const T* x, const int* kk, T* yy) { t_slope(tab, x, kk, yy); };
    auto cl = [&](const T* x, T* yy) { c_left(*cache, x, yy); };
    auto cn = [&](const T* x, T* yy) { c_naive(*cache, x, yy); };
    run_table(tl); same("table L");
    run_cached(cl); same("cached L");
    std::printf("lanes changing interval per chunk: %.2f %%\n", 100.0 * moved / P);
#if defined(__FMA__)
    auto tf = [&](const T* x, const int* kk, T* yy) { t_left_fma(tab, x, kk, yy); };
    auto cf = [&](const T* x, T* yy) { c_left_fma(*cache, x, yy); };
    run_table(tf); same("table L fma");
    run_cached(cf); same("cached L fma");
#endif

    auto nop = [](const T*, T*) {};
    const double t_naive = best_ns([&] { run_table(ts); }, P);
    const double t_l = best_ns([&] { run_table(tl); }, P);
    const double c_nv = best_ns([&] { run_cached(cn); }, P);
    const double c_l = best_ns([&] { run_cached(cl); }, P);
    const double c_ref = best_ns([&] { run_cached(nop); }, P);
    double t_f = -1, c_f = -1;
#if defined(__FMA__)
    t_f = best_ns([&] { run_table(tf); }, P);
    c_f = best_ns([&] { run_cached(cf); }, P);
#endif
    auto f = [](double v) { return v > 0 ? report::num(v, "%.2f") : std::string("-"); };
    report::Table t({"layout", "naive slope", "L", "L fma", "L / naive"}, 12);
    t.line({"table", f(t_naive), f(t_l), f(t_f), report::num((t_f > 0 ? t_f : t_l) / t_naive, "%.2fx")});
    t.line({"cached", f(c_nv), f(c_l), f(c_f), report::num((c_f > 0 ? c_f : c_l) / c_nv, "%.2fx")});
    std::printf("(ns per point; cached refresh alone: %.2f ns; L / naive uses L fma when available)\n", c_ref);
    report::Csv csv(std::string("t12_chunk_opt_") + fp::type_name<T>() + "_" + label,
                    "layout,naive_slope_ns,L_ns,L_fma_ns,refresh_ns");
    csv.row("table", t_naive, t_l, t_f, 0);
    csv.row("cached", c_nv, c_l, c_f, c_ref);
    delete cache;
}

int main() {
#if defined(__FMA__)
    std::printf("t12: optimizing L in chunks of %d lanes (FMA available)\n", M);
#else
    std::printf("t12: optimizing L in chunks of %d lanes (no FMA in this build)\n", M);
#endif
    report::Checks checks;
    for (auto [step, label] : {std::pair{0.002, "slow"}, {0.05, "moving"}}) {
        run<float>(checks, step, label);
        run<double>(checks, step, label);
    }
    return checks.exit_code();
}
