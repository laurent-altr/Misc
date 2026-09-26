// Property 6: behaviour at the interior knots of a table. The lookup puts
// x_k in interval [x_k, x_k+1], so the value just below the knot comes from
// the previous interval. Measured on many random tables:
//   left!=Yk : f_{k-1}(x_k) != Y[k]   (left interval does not reach the knot)
//   right!=Yk: f_k(x_k) != Y[k]
//   jump     : max |f_{k-1}(x_k) - f_k(x_k)| in ulps of max(|Y[k-1]|,|Y[k]|,|Y[k+1]|)
//   mono     : on increasing tables, violations of
//              f(prev(x_k)) <= f(x_k) <= f(next(x_k)) through the lookup
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "cases.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "report.hpp"

template <class T> T lookup(interp::Fn<T> f, const cases::Table<T>& t, T x) {
    auto it = std::upper_bound(t.x.begin(), t.x.end(), x);
    std::size_t k = std::clamp<std::ptrdiff_t>(it - t.x.begin() - 1, 0, t.x.size() - 2);
    return f(t.x[k], t.x[k + 1], t.y[k], t.y[k + 1], x);
}

template <class T> void run(report::Checks& checks) {
    std::printf("\n=== %s ===\n", fp::type_name<T>());
    auto any = cases::random_tables<T>(2000, 40, false);
    auto inc = cases::random_tables<T>(2000, 40, true, 999);

    report::Table tab({"formula", "knots", "left!=Yk", "right!=Yk", "max jump", "mono viol"}, 10);
    report::Csv csv(std::string("t05_knots_") + fp::type_name<T>(), "formula,knots,left_ne,right_ne,max_jump_ulps,mono_checks,mono_violations");
    for (const auto& fm : interp::formulas<T>()) {
        long knots = 0, left_ne = 0, right_ne = 0, mono_n = 0, mono_bad = 0;
        double jump = 0;
        for (const auto& tabs : {&any, &inc}) {
            for (const auto& t : *tabs) {
                for (std::size_t k = 1; k + 1 < t.x.size(); ++k) {
                    T l = fm.f(t.x[k - 1], t.x[k], t.y[k - 1], t.y[k], t.x[k]);
                    T r = fm.f(t.x[k], t.x[k + 1], t.y[k], t.y[k + 1], t.x[k]);
                    ++knots;
                    left_ne += l != t.y[k];
                    right_ne += r != t.y[k];
                    T scale = std::max({std::fabs(t.y[k - 1]), std::fabs(t.y[k]), std::fabs(t.y[k + 1])});
                    jump = std::max(jump, std::fabs(double(l) - double(r)) / double(fp::ulp(scale)));
                    if (tabs == &inc) {
                        T a = lookup(fm.f, t, fp::next_down(t.x[k]));
                        T b = lookup(fm.f, t, t.x[k]);
                        T c = lookup(fm.f, t, fp::next_up(t.x[k]));
                        ++mono_n;
                        mono_bad += !(a <= b && b <= c);
                    }
                }
            }
        }
        tab.line({fm.id, std::to_string(knots), std::to_string(left_ne), std::to_string(right_ne), report::num(jump),
                  std::to_string(mono_bad) + "/" + std::to_string(mono_n)});
        csv.row(fm.name, knots, left_ne, right_ne, jump, mono_n, mono_bad);
        if (std::string(fm.id) == "F")
            checks.expect(left_ne == 0 && right_ne == 0 && mono_bad == 0,
                          std::string("F (") + fp::type_name<T>() + ") continuous and monotonic at knots");
    }
}

int main() {
    std::printf("t05: continuity and monotonicity at the knots of random tables\n");
    report::Checks checks;
    run<float>(checks);
    run<double>(checks);
    return checks.exit_code();
}
