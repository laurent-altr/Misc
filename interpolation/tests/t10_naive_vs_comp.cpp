// Head-to-head: the naive formula anchored on x0 (A) against the fma-free
// compensated nearest-bound formula (H'), on numerical quality and speed.
// Speed is measured with the formulas inlined (no function pointer), in three
// usage patterns, plus H' with per-segment data precomputed (Hp).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>

#include "cases.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "reference.hpp"
#include "report.hpp"

struct Naive {
    template <class T> T operator()(T x0, T x1, T y0, T y1, T x) const { return interp::left(x0, x1, y0, y1, x); }
};
struct Comp {
    template <class T> T operator()(T x0, T x1, T y0, T y1, T x) const {
        return interp::nearest_comp_nofma(x0, x1, y0, y1, x);
    }
};

struct Quality {
    long n = 0, cr = 0, steps = 0, same_step = 0, stalls = 0, backwards = 0, segs = 0, ex0 = 0, ex1 = 0, knots = 0, knot_ne = 0;
    std::int64_t max_ulps = 0, extra_flat = 0;
    double sum_ulps = 0;
};

template <class T, class F> Quality quality(F f) {
    Quality q;
    auto segs = cases::random<T>(2000, 31);
    std::mt19937_64 g(5);
    std::uniform_real_distribution<double> u(0, 1);
    auto flat_run = [](const std::vector<T>& y) {
        std::int64_t run = 0, best = 0;
        for (std::size_t i = 1; i < y.size(); ++i) best = std::max(best, run = y[i] == y[i - 1] ? run + 1 : 0);
        return best;
    };
    for (const auto& s : segs) {
        ++q.segs;
        q.ex0 += f(s.x0, s.x1, s.y0, s.y1, s.x0) == s.y0;
        q.ex1 += f(s.x0, s.x1, s.y0, s.y1, s.x1) == s.y1;
        // Pointwise, random x.
        for (int i = 0; i < 512; ++i) {
            T x = std::clamp(T(s.x0 + u(g) * (s.x1 - s.x0)), s.x0, s.x1);
            T y = f(s.x0, s.x1, s.y0, s.y1, x), c = ref::rounded(s.x0, s.x1, s.y0, s.y1, x);
            std::int64_t d = std::abs(fp::ulp_distance(c, y));
            ++q.n;
            q.cr += d == 0;
            q.sum_ulps += double(d);
            q.max_ulps = std::max(q.max_ulps, d);
        }
        // Runs of 512 consecutive x at a random t and near x1.
        const bool inc = s.y1 >= s.y0;
        for (double t0 : {u(g), 0.97}) {
            std::vector<T> y, c;
            for (T x = T(s.x0 + t0 * (s.x1 - s.x0)); x < s.x1 && y.size() < 513; x = fp::next_up(x)) {
                y.push_back(f(s.x0, s.x1, s.y0, s.y1, x));
                c.push_back(ref::rounded(s.x0, s.x1, s.y0, s.y1, x));
            }
            for (std::size_t i = 1; i < y.size(); ++i) {
                std::int64_t df = fp::ulp_distance(y[i - 1], y[i]), dc = fp::ulp_distance(c[i - 1], c[i]);
                ++q.steps;
                q.same_step += df == dc;
                q.stalls += df == 0 && dc != 0;
                q.backwards += inc ? df < 0 : df > 0;
            }
            q.extra_flat = std::max(q.extra_flat, flat_run(y) - flat_run(c));
        }
    }
    // Knots of random tables: does the interval on the left reach Y[k]?
    for (const auto& t : cases::random_tables<T>(1000, 40, false)) {
        for (std::size_t k = 1; k + 1 < t.x.size(); ++k) {
            ++q.knots;
            q.knot_ne += f(t.x[k - 1], t.x[k], t.y[k - 1], t.y[k], t.x[k]) != t.y[k];
        }
    }
    return q;
}

template <class T> void print_quality() {
    Quality a = quality<T>(Naive{}), h = quality<T>(Comp{});
    report::Table tab({"metric", "A (naive)", "H'"}, 12);
    auto pct = [](long x, long n) { return report::num(100.0 * x / n, "%.2f%%"); };
    report::Csv csv(std::string("t10_quality_") + fp::type_name<T>(), "metric,naive,comp_nofma");
    auto line = [&](const char* m, const std::string& va, const std::string& vh) {
        tab.line({m, va, vh});
        csv.row(m, va, vh);
    };
    line("correctly rounded", pct(a.cr, a.n), pct(h.cr, h.n));
    line("mean err (ulp)", report::num(a.sum_ulps / a.n), report::num(h.sum_ulps / h.n));
    line("max err (ulp)", report::num(double(a.max_ulps)), report::num(double(h.max_ulps)));
    line("steps = CR", pct(a.same_step, a.steps), pct(h.same_step, h.steps));
    line("stalls", std::to_string(a.stalls), std::to_string(h.stalls));
    line("backward steps", std::to_string(a.backwards), std::to_string(h.backwards));
    line("extra flat run", std::to_string(a.extra_flat), std::to_string(h.extra_flat));
    line("exact at x0", pct(a.ex0, a.segs), pct(h.ex0, h.segs));
    line("exact at x1", pct(a.ex1, a.segs), pct(h.ex1, h.segs));
    line("knot mismatch", pct(a.knot_ne, a.knots), pct(h.knot_ne, h.knots));
}

// ---------------------------------------------------------------- speed
// Each pattern is an out-of-line loop over plain arrays (restrict pointers),
// so the compiler may inline and vectorize the formula as in real code.

template <class Fn> double best_ns(Fn&& fn, long calls, int reps = 100) {
    double best = 1e30;
    for (int r = 0; r < reps; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        fn();
        auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count() / calls);
    }
    return best;
}

template <class T> struct In { T x0, x1, y0, y1, x; };

// 1. Independent calls: a different segment and x each time.
template <class T, class F>
__attribute__((noinline)) void loop_indep(const In<T>* __restrict e, T* __restrict y, int n, F f) {
    for (int i = 0; i < n; ++i) y[i] = f(e[i].x0, e[i].x1, e[i].y0, e[i].y1, e[i].x);
}

// 2. Table: binary search, then interpolation.
template <class T> std::size_t find(const T* X, int n, T x) {
    auto k = std::upper_bound(X, X + n, x) - X - 1;
    return std::size_t(std::clamp<std::ptrdiff_t>(k, 0, n - 2));
}
template <class T, class F>
__attribute__((noinline)) void loop_table(const T* __restrict X, const T* __restrict Y, int nk, const T* __restrict q,
                                          T* __restrict y, int n, F f) {
    for (int i = 0; i < n; ++i) {
        std::size_t k = find(X, nk, q[i]);
        y[i] = f(X[k], X[k + 1], Y[k], Y[k + 1], q[i]);
    }
}
template <class T>
__attribute__((noinline)) void loop_table_pre(const T* __restrict X, const interp::CompSeg<T>* __restrict pre, int nk,
                                              const T* __restrict q, T* __restrict y, int n) {
    for (int i = 0; i < n; ++i) y[i] = interp::eval_comp_seg(pre[find(X, nk, q[i])], q[i]);
}
template <class T>
__attribute__((noinline)) void loop_search(const T* __restrict X, const T* __restrict Y, int nk, const T* __restrict q,
                                           T* __restrict y, int n) {
    for (int i = 0; i < n; ++i) y[i] = Y[find(X, nk, q[i])];
}

// 2b. Hunt: queries move a little each time; the interval is found by walking
// from the previous one (cheap, as with a good initial guess).
template <class T> std::size_t hunt(const T* X, int n, T x, std::size_t k) {
    while (k > 0 && x < X[k]) --k;
    while (k + 2 < std::size_t(n) && x >= X[k + 1]) ++k;
    return k;
}
template <class T, class F>
__attribute__((noinline)) void loop_hunt(const T* __restrict X, const T* __restrict Y, int nk, const T* __restrict q,
                                         T* __restrict y, int n, F f) {
    std::size_t k = 0;
    for (int i = 0; i < n; ++i) {
        k = hunt(X, nk, q[i], k);
        y[i] = f(X[k], X[k + 1], Y[k], Y[k + 1], q[i]);
    }
}
template <class T>
__attribute__((noinline)) void loop_hunt_pre(const T* __restrict X, const interp::CompSeg<T>* __restrict pre, int nk,
                                             const T* __restrict q, T* __restrict y, int n) {
    std::size_t k = 0;
    for (int i = 0; i < n; ++i) {
        k = hunt(X, nk, q[i], k);
        y[i] = interp::eval_comp_seg(pre[k], q[i]);
    }
}
template <class T>
__attribute__((noinline)) void loop_hunt_only(const T* __restrict X, const T* __restrict Y, int nk,
                                              const T* __restrict q, T* __restrict y, int n) {
    std::size_t k = 0;
    for (int i = 0; i < n; ++i) y[i] = Y[k = hunt(X, nk, q[i], k)];
}

// 3. Sweep: one segment, an array of x -> an array of y.
template <class T, class F>
__attribute__((noinline)) void loop_sweep(T a0, T a1, T b0, T b1, const T* __restrict x, T* __restrict y, int n, F f) {
    for (int i = 0; i < n; ++i) y[i] = f(a0, a1, b0, b1, x[i]);
}
template <class T>
__attribute__((noinline)) void loop_sweep_pre(interp::CompSeg<T> g, const T* __restrict x, T* __restrict y, int n) {
    for (int i = 0; i < n; ++i) y[i] = interp::eval_comp_seg(g, x[i]);
}

template <class T> void print_speed() {
    std::mt19937_64 g(9);
    std::uniform_real_distribution<double> u(0, 1), v(-100, 100);
    constexpr int N = 1 << 14;
    std::vector<T> out(N);

    std::vector<In<T>> in(N);
    for (auto& e : in) {
        e.x0 = T(v(g)); e.x1 = T(e.x0 + 0.001 + 10 * u(g)); e.y0 = T(v(g)); e.y1 = T(v(g));
        e.x = T(e.x0 + u(g) * (e.x1 - e.x0));
    }
    auto indep = [&](auto f) { return best_ns([&] { loop_indep(in.data(), out.data(), N, f); }, N); };

    std::vector<T> X, Y;
    { double x = 0; for (int k = 0; k < 1024; ++k) { X.push_back(T(x)); Y.push_back(T(v(g))); x += 0.01 + u(g); } }
    const int nk = int(X.size());
    std::vector<interp::CompSeg<T>> pre;
    for (int k = 0; k + 1 < nk; ++k) pre.push_back(interp::make_comp_seg(X[k], X[k + 1], Y[k], Y[k + 1]));
    std::vector<T> q(N);
    for (auto& x : q) x = T(u(g) * X.back());
    auto table = [&](auto f) { return best_ns([&] { loop_table(X.data(), Y.data(), nk, q.data(), out.data(), N, f); }, N); };
    const double table_pre = best_ns([&] { loop_table_pre(X.data(), pre.data(), nk, q.data(), out.data(), N); }, N);
    const double search = best_ns([&] { loop_search(X.data(), Y.data(), nk, q.data(), out.data(), N); }, N);

    // Random walk over the table: steps of up to 5 % of the mean knot spacing,
    // so about one query in 20 moves to a neighbouring interval.
    std::vector<T> w(N);
    {
        const double span = double(X.back()), mean = span / (nk - 1);
        double x = span / 2;
        for (auto& e : w) {
            x = std::clamp(x + (u(g) - 0.5) * 0.1 * mean, 0.0, span * 0.999);
            e = T(x);
        }
    }
    auto walk = [&](auto f) { return best_ns([&] { loop_hunt(X.data(), Y.data(), nk, w.data(), out.data(), N, f); }, N); };
    const double walk_pre = best_ns([&] { loop_hunt_pre(X.data(), pre.data(), nk, w.data(), out.data(), N); }, N);
    const double walk_only = best_ns([&] { loop_hunt_only(X.data(), Y.data(), nk, w.data(), out.data(), N); }, N);

    std::vector<T> xs(N);
    const T a0 = T(1.3), a1 = T(7.9), b0 = T(-4.2), b1 = T(11.7);
    for (int i = 0; i < N; ++i) xs[i] = T(a0 + (a1 - a0) * i / N);
    auto sweep = [&](auto f) { return best_ns([&] { loop_sweep(a0, a1, b0, b1, xs.data(), out.data(), N, f); }, N); };
    const auto g1 = interp::make_comp_seg(a0, a1, b0, b1);
    const double sweep_pre = best_ns([&] { loop_sweep_pre(g1, xs.data(), out.data(), N); }, N);

    report::Table tab({"pattern", "A (naive)", "H'", "Hp (pre)", "H'/A", "Hp/A"}, 10);
    report::Csv csv(std::string("t10_speed_") + fp::type_name<T>(), "pattern,naive_ns,comp_ns,comp_pre_ns");
    auto line = [&](const char* name, double a, double h, double hp) {
        tab.line({name, report::num(a, "%.2f"), report::num(h, "%.2f"), hp > 0 ? report::num(hp, "%.2f") : "-",
                  report::num(h / a, "%.1fx"), hp > 0 ? report::num(hp / a, "%.1fx") : "-"});
        csv.row(name, a, h, hp);
    };
    line("independent", indep(Naive{}), indep(Comp{}), -1);
    line("table", table(Naive{}), table(Comp{}), table_pre);
    line("hunt", walk(Naive{}), walk(Comp{}), walk_pre);
    line("sweep", sweep(Naive{}), sweep(Comp{}), sweep_pre);
    std::printf("(ns per call; the search alone costs %.2f ns in the table pattern, %.2f ns in the hunt pattern)\n",
                search, walk_only);
}

template <class T> void check_pre(report::Checks& checks) {
    long diff = 0;
    for (const auto& s : cases::random<T>(500, 8)) {
        const auto g = interp::make_comp_seg(s.x0, s.x1, s.y0, s.y1);
        for (T x = s.x0 + (s.x1 - s.x0) / 3, i = 0; i < 200; ++i, x = fp::next_up(x))
            diff += fp::ordered(interp::eval_comp_seg(g, x)) != fp::ordered(interp::nearest_comp_nofma(s.x0, s.x1, s.y0, s.y1, x));
    }
    checks.expect(diff == 0, std::string("precomputed H' identical to H' (") + fp::type_name<T>() + ")");
}

int main() {
    std::printf("t10: naive formula A vs compensated H' (fma-free)\n");
    report::Checks checks;
    std::printf("\n=== quality, float ===\n");
    print_quality<float>();
    std::printf("\n=== quality, double ===\n");
    print_quality<double>();
    std::printf("\n=== speed, float ===\n");
    print_speed<float>();
    std::printf("\n=== speed, double ===\n");
    print_speed<double>();
    check_pre<float>(checks);
    check_pre<double>(checks);
    return checks.exit_code();
}
