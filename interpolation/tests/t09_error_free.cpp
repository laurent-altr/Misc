// Do the error-free transformations survive the compiler? two_sum and
// two_prod must satisfy a + b == s + e and a * b == p + e exactly (checked in
// __float128). Under -ffast-math the compiler may simplify e to 0; under FMA
// contraction it may fuse the multiply inside Veltkamp's split. Also compares
// H (fma) and H' (no fma) bit for bit.
#include <cmath>
#include <cstdio>
#include <random>

#include "fp_utils.hpp"
#include "interp.hpp"
#include "reference.hpp"
#include "report.hpp"

template <class T> void run(report::Checks& checks) {
    std::printf("\n=== %s ===\n", fp::type_name<T>());
    std::mt19937_64 g(7);
    std::uniform_real_distribution<double> m(-1, 1), ex(-20, 20);
    auto rnd = [&] { return T(std::ldexp(m(g), int(ex(g)))); };
    constexpr long N = 1 << 21;
    long sum_bad = 0, prod_bad = 0, sum_zero = 0, prod_zero = 0, h_diff = 0;
    for (long i = 0; i < N; ++i) {
        T a = rnd(), b = rnd(), s, e;
        interp::two_sum(a, b, s, e);
        sum_bad += ref::Q(s) + ref::Q(e) != ref::Q(a) + ref::Q(b);
        sum_zero += e == 0;
        interp::two_prod(a, b, s, e);
        prod_bad += ref::Q(s) + ref::Q(e) != ref::Q(a) * ref::Q(b);
        prod_zero += e == 0;

        T x0 = T(100 * m(g)), x1 = T(x0 + 0.001 + 10 * std::fabs(m(g))), y0 = T(100 * m(g)), y1 = T(100 * m(g));
        T x = T(x0 + std::fabs(m(g)) * (x1 - x0));
        T h1 = interp::nearest_comp(x0, x1, y0, y1, x), h2 = interp::nearest_comp_nofma(x0, x1, y0, y1, x);
        h_diff += fp::ordered(h1) != fp::ordered(h2);
    }
    report::Table tab({"check", "failures", "e == 0"}, 12);
    tab.line({"two_sum", std::to_string(sum_bad), std::to_string(sum_zero)});
    tab.line({"two_prod", std::to_string(prod_bad), std::to_string(prod_zero)});
    tab.line({"H vs H'", std::to_string(h_diff), "-"});
    std::printf("(%ld random cases each)\n", N);
    report::Csv csv(std::string("t09_error_free_") + fp::type_name<T>(), "check,cases,failures,zero_error_terms");
    csv.row("two_sum", N, sum_bad, sum_zero);
    csv.row("two_prod", N, prod_bad, prod_zero);
    csv.row("h_vs_hnofma", N, h_diff, 0);
    const std::string t = fp::type_name<T>();
    checks.expect(sum_bad == 0, "two_sum exact (" + t + ")");
    checks.expect(prod_bad == 0, "two_prod exact (" + t + ")");
}

int main() {
    std::printf("t09: error-free transformations under the current compiler flags\n");
    report::Checks checks;
    run<float>(checks);
    run<double>(checks);
    return checks.exit_code();
}
