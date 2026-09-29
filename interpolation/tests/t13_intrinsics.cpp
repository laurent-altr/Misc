// Hand-written AVX2 / AVX-512 intrinsics (immintrin.h) against the
// compiler-vectorized cached L with FMA (chunk::eval_cached_fma), double,
// chunks of 128 lanes. Same arithmetic, same bits (checked). Also a refresh
// using a comparison bit mask to visit only the lanes that left their interval.
// Built only when the target has AVX2 + FMA (AVX-512 parts need AVX512F).
#if defined(__AVX2__) && defined(__FMA__)
#include <immintrin.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include "chunk.hpp"
#include "fp_utils.hpp"
#include "report.hpp"
constexpr int M = 128, C = 512;
using Ln = chunk::Lanes<double, M>;

__attribute__((noinline)) void eval_compiler(const Ln& c, const double* x, double* y) { chunk::eval_cached_fma(c, x, y); }

// Same arithmetic, AVX2 intrinsics (4 doubles).
__attribute__((noinline)) void eval_avx2(const Ln& c, const double* __restrict x, double* __restrict y) {
    for (int i = 0; i < M; i += 4) {
        const __m256d xi = _mm256_loadu_pd(x + i), x0 = _mm256_load_pd(c.x0 + i), y0 = _mm256_load_pd(c.y0 + i);
        const __m256d s = _mm256_load_pd(c.s + i), sl = _mm256_load_pd(c.s_lo + i);
        const __m256d nx0 = _mm256_sub_pd(_mm256_setzero_pd(), x0);
        // two_sum(x, -x0)
        const __m256d h = _mm256_add_pd(xi, nx0);
        const __m256d bb = _mm256_sub_pd(h, xi);
        const __m256d he = _mm256_add_pd(_mm256_sub_pd(xi, _mm256_sub_pd(h, bb)), _mm256_sub_pd(nx0, bb));
        const __m256d p = _mm256_mul_pd(h, s);
        const __m256d pe = _mm256_fmsub_pd(h, s, p);
        // two_sum(y0, p)
        const __m256d sum = _mm256_add_pd(y0, p);
        const __m256d b2 = _mm256_sub_pd(sum, y0);
        const __m256d se = _mm256_add_pd(_mm256_sub_pd(y0, _mm256_sub_pd(sum, b2)), _mm256_sub_pd(p, b2));
        const __m256d corr = _mm256_add_pd(_mm256_mul_pd(h, sl), _mm256_mul_pd(he, s));
        _mm256_storeu_pd(y + i, _mm256_add_pd(sum, _mm256_add_pd(se, _mm256_add_pd(pe, corr))));
    }
}

#if defined(__AVX512F__)
// AVX-512 (8 doubles).
__attribute__((noinline)) void eval_avx512(const Ln& c, const double* __restrict x, double* __restrict y) {
    for (int i = 0; i < M; i += 8) {
        const __m512d xi = _mm512_loadu_pd(x + i), x0 = _mm512_load_pd(c.x0 + i), y0 = _mm512_load_pd(c.y0 + i);
        const __m512d s = _mm512_load_pd(c.s + i), sl = _mm512_load_pd(c.s_lo + i);
        const __m512d nx0 = _mm512_sub_pd(_mm512_setzero_pd(), x0);
        const __m512d h = _mm512_add_pd(xi, nx0);
        const __m512d bb = _mm512_sub_pd(h, xi);
        const __m512d he = _mm512_add_pd(_mm512_sub_pd(xi, _mm512_sub_pd(h, bb)), _mm512_sub_pd(nx0, bb));
        const __m512d p = _mm512_mul_pd(h, s);
        const __m512d pe = _mm512_fmsub_pd(h, s, p);
        const __m512d sum = _mm512_add_pd(y0, p);
        const __m512d b2 = _mm512_sub_pd(sum, y0);
        const __m512d se = _mm512_add_pd(_mm512_sub_pd(y0, _mm512_sub_pd(sum, b2)), _mm512_sub_pd(p, b2));
        const __m512d corr = _mm512_add_pd(_mm512_mul_pd(h, sl), _mm512_mul_pd(he, s));
        _mm512_storeu_pd(y + i, _mm512_add_pd(sum, _mm512_add_pd(se, _mm512_add_pd(pe, corr))));
    }
}

#else
#define eval_avx512 eval_avx2
#endif

// Refresh check with intrinsics: movemask, early exit.
__attribute__((noinline)) int refresh_avx2(const chunk::Table<double>& t, Ln& c, const double* x) {
    __m256d any = _mm256_setzero_pd();
    for (int i = 0; i < M; i += 4) {
        const __m256d xi = _mm256_loadu_pd(x + i);
        any = _mm256_or_pd(any, _mm256_or_pd(_mm256_cmp_pd(xi, _mm256_load_pd(c.lo + i), _CMP_LT_OQ),
                                             _mm256_cmp_pd(xi, _mm256_load_pd(c.hi + i), _CMP_GE_OQ)));
    }
    if (!_mm256_movemask_pd(any)) return 0;
    int moved = 0;
    for (int i = 0; i < M; i += 4) {  // visit only the flagged lanes via the bit mask
        const __m256d xi = _mm256_loadu_pd(x + i);
        int m = _mm256_movemask_pd(_mm256_or_pd(_mm256_cmp_pd(xi, _mm256_load_pd(c.lo + i), _CMP_LT_OQ),
                                                _mm256_cmp_pd(xi, _mm256_load_pd(c.hi + i), _CMP_GE_OQ)));
        while (m) {
            const int l = i + __builtin_ctz(m);
            chunk::fill_lane(t, c, l, chunk::hunt(t, x[l], c.k[l]));
            ++moved;
            m &= m - 1;
        }
    }
    return moved;
}
__attribute__((noinline)) int refresh_compiler(const chunk::Table<double>& t, Ln& c, const double* x) { return chunk::refresh(t, c, x); }

template <class F> double best(F f, long n) {
    double b = 1e30;
    for (int r = 0; r < 100; ++r) {
        auto t0 = std::chrono::steady_clock::now(); f();
        b = std::min(b, std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / n);
    }
    return b;
}

int main() {
    std::printf("t13: intrinsics vs compiler vectorization (double, %d lanes, %s)\n", M,
#if defined(__AVX512F__)
                "AVX-512 available"
#else
                "no AVX-512: the AVX-512 column repeats AVX2"
#endif
    );
    report::Checks checks;
    report::Csv csv("t13_intrinsics", "step,compiler_ns,avx2_ns,avx512_ns,refresh_compiler_ns,refresh_intr_ns");
    std::mt19937_64 g(21); std::uniform_real_distribution<double> u(0, 1), v(-100, 100);
    std::vector<double> X, Y; double xx = 0;
    for (int k = 0; k < 1024; ++k) { X.push_back(xx); Y.push_back(v(g)); xx += 0.01 + u(g); }
    auto tab = chunk::make_table(X, Y);
    for (double step : {0.002, 0.05}) {
        const double span = X.back(), mean = span / 1023;
        std::vector<double> xs(size_t(C) * M); std::vector<int> k0(M);
        for (int i = 0; i < M; ++i) {
            double x = u(g) * span * 0.999;
            k0[i] = int(std::upper_bound(X.begin(), X.end(), x) - X.begin()) - 1;
            for (int c = 0; c < C; ++c) { x = std::clamp(x + (u(g) - 0.5) * 2 * step * mean, 0.0, span * 0.999); xs[size_t(c) * M + i] = x; }
        }
        std::vector<double> y1(xs.size()), y2(xs.size()), y3(xs.size());
        auto* cache = new Ln;
        auto run = [&](auto refresh, auto eval, std::vector<double>& y) {
            chunk::init(tab, *cache, &xs[0], k0.data());
            for (int c = 0; c < C; ++c) { refresh(tab, *cache, &xs[size_t(c) * M]); eval(*cache, &xs[size_t(c) * M], &y[size_t(c) * M]); }
        };
        run(refresh_compiler, eval_compiler, y1); run(refresh_avx2, eval_avx2, y2); run(refresh_avx2, eval_avx512, y3);
        long d = 0; for (size_t i = 0; i < y1.size(); ++i) d += fp::ordered(y1[i]) != fp::ordered(y2[i]) || fp::ordered(y1[i]) != fp::ordered(y3[i]);
        checks.expect(d == 0, "intrinsics give the same bits as the compiler version");
        const long P = long(C) * M;
        auto nop = [](const Ln&, const double*, double*) {};
        const double a = best([&] { run(refresh_compiler, eval_compiler, y1); }, P);
        const double b = best([&] { run(refresh_avx2, eval_avx2, y2); }, P);
        const double e = best([&] { run(refresh_avx2, eval_avx512, y3); }, P);
        const double rc = best([&] { run(refresh_compiler, nop, y1); }, P);
        const double ri = best([&] { run(refresh_avx2, nop, y1); }, P);
        std::printf("\nx steps up to %.1f %% of the knot spacing (refresh + evaluation, ns per point)\n", 100 * step);
        report::Table t({"variant", "total", "refresh only"}, 12);
        t.line({"compiler", report::num(a, "%.2f"), report::num(rc, "%.2f")});
        t.line({"AVX2 intr.", report::num(b, "%.2f"), report::num(ri, "%.2f")});
        t.line({"AVX-512 intr.", report::num(e, "%.2f"), report::num(ri, "%.2f")});
        csv.row(step, a, b, e, rc, ri);
        delete cache;
    }
    return checks.exit_code();
}
#else
#include <cstdio>
int main() {
    std::printf("t13: skipped (build without AVX2 + FMA, e.g. add -march=native)\n");
    return 0;
}
#endif
