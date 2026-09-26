// Property 3: monotonicity when x moves by one ulp. For an increasing
// segment, f(next_up(x)) >= f(x) must hold (<= for a decreasing one).
// Formulas made only of correctly rounded operations that are each monotone
// in x (A, A', A'', B, D, F) cannot violate it; C (sum of an increasing and a
// decreasing term) and E (switches formula at mid-interval) can.
#include <cstdio>
#include <map>

#include "cases.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "report.hpp"
#include "sampling.hpp"

template <class T> void run(report::Checks& checks) {
    std::printf("\n=== %s ===  violations (max size in ulps of y)\n", fp::type_name<T>());
    auto segs = cases::fixed<T>();
    auto rnd = cases::random<T>(16);
    for (auto& s : rnd) s.name = "rand";
    segs.insert(segs.end(), rnd.begin(), rnd.end());

    std::vector<std::string> names;  // column order, "rand" aggregated
    for (const auto& s : segs)
        if (names.empty() || names.back() != s.name) names.push_back(s.name);

    auto fms = interp::formulas<T>();
    struct Stat { long viol = 0, n = 0; std::int64_t worst = 0; };
    std::map<std::string, std::map<std::string, Stat>> st;  // formula -> case -> stat

    sampling::Budget budget;
    for (const auto& s : segs) {
        const bool inc = s.y1 >= s.y0;
        sampling::steps(s, budget, [&](T x) {
            T xp = fp::next_up(x);
            for (const auto& fm : fms) {
                T a = fm.f(s.x0, s.x1, s.y0, s.y1, x), b = fm.f(s.x0, s.x1, s.y0, s.y1, xp);
                auto& c = st[fm.id][s.name];
                ++c.n;
                std::int64_t back = inc ? fp::ulp_distance(b, a) : fp::ulp_distance(a, b);
                if (back > 0) {
                    ++c.viol;
                    c.worst = std::max(c.worst, back);
                }
            }
        });
    }

    std::vector<std::string> head{"formula"};
    head.insert(head.end(), names.begin(), names.end());
    report::Table tab(head, 10);
    report::Csv csv(std::string("t02_monotonic_") + fp::type_name<T>(), "formula,case,steps,violations,worst_ulps");
    for (const auto& fm : fms) {
        std::vector<std::string> row{fm.id};
        long total = 0;
        for (const auto& n : names) {
            const Stat& c = st[fm.id][n];
            row.push_back(c.viol ? std::to_string(c.viol) + " (" + std::to_string(c.worst) + ")" : "0");
            csv.row(fm.name, n, c.n, c.viol, c.worst);
            total += c.viol;
        }
        tab.line(row);
        const std::string id = fm.id;
        if (id != "C" && id != "E")
            checks.expect(total == 0, id + " (" + fp::type_name<T>() + ") monotonic");
    }
    long steps = 0;
    for (const auto& n : names) steps += st["A"][n].n;
    std::printf("(%ld steps per formula)\n", steps);
}

int main() {
    std::printf("t02: monotonicity under a one-ulp move of x\n");
    report::Checks checks;
    run<float>(checks);
    run<double>(checks);
    return checks.exit_code();
}
