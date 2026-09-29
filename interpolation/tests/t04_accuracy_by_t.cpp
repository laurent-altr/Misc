// Property 5: pointwise error of f(x) against the exact value, by position
// t = (x - x0) / (x1 - x0) in ten bins. Errors are in ulps of the segment
// scale max(|y0|,|y1|); the last column gives the maximum in ulps of y itself
// (large near a zero crossing, where the relative error explodes).
#include <array>
#include <cmath>
#include <cstdio>

#include "cases.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "reference.hpp"
#include "report.hpp"
#include "sampling.hpp"

template <class T> void run(report::Checks& checks) {
    std::printf("\n=============== %s ===============\n", fp::type_name<T>());
    auto segs = cases::fixed<T>();
    auto rnd = cases::random<T>(64);
    for (auto& s : rnd) s.name = "rand";
    segs.insert(segs.end(), rnd.begin(), rnd.end());

    auto fms = interp::formulas<T>();
    const std::size_t nf = fms.size() + 1;  // + CR baseline
    constexpr int B = 10;
    report::Csv csv(std::string("t04_accuracy_by_t_") + fp::type_name<T>(), "case,formula,bin_lo,bin_hi,samples,max_err_ulps,mean_err_ulps");

    std::size_t i = 0;
    while (i < segs.size()) {
        const std::string name = segs[i].name;
        struct Bin { long n = 0; double max = 0, sum = 0; };
        std::vector<std::array<Bin, B>> st(nf);
        std::vector<double> local(nf, 0);
        for (; i < segs.size() && segs[i].name == name; ++i) {
            const auto& s = segs[i];
            const ref::Q u = fp::ulp(std::max(std::fabs(s.y0), std::fabs(s.y1)));
            sampling::uniform(s, std::int64_t(1) << 18, [&](T x) {
                ref::Q e = ref::value(s.x0, s.x1, s.y0, s.y1, x);
                double t = static_cast<double>((ref::Q(x) - s.x0) / (ref::Q(s.x1) - s.x0));
                int b = std::min(B - 1, static_cast<int>(t * B));
                for (std::size_t k = 0; k < nf; ++k) {
                    T y = k < fms.size() ? fms[k].f(s.x0, s.x1, s.y0, s.y1, x) : static_cast<T>(e);
                    double err = std::fabs(static_cast<double>((ref::Q(y) - e) / u));
                    Bin& bin = st[k][b];
                    ++bin.n;
                    bin.sum += err;
                    bin.max = std::max(bin.max, err);
                    local[k] = std::max(local[k], std::fabs(ref::err_ulps(y, e)));
                }
            });
        }

        std::printf("\n-- case %s: max error per t bin, ulps of max(|y0|,|y1|)\n", name.c_str());
        std::vector<std::string> head{"formula"};
        for (int b = 0; b < B; ++b) head.push_back(report::num(b / double(B), "%.1f") + "-");
        head.push_back("ulp(y)");
        report::Table tab(head, 7);
        for (std::size_t k = 0; k < nf; ++k) {
            const std::string id = k < fms.size() ? fms[k].id : "CR";
            std::vector<std::string> row{id};
            for (int b = 0; b < B; ++b) {
                const Bin& bin = st[k][b];
                row.push_back(bin.n ? report::num(bin.max, "%.2f") : "-");
                csv.row(name, k < fms.size() ? fms[k].name : "correctly_rounded", b / double(B), (b + 1) / double(B),
                        bin.n, bin.max, bin.n ? bin.sum / bin.n : 0);
            }
            row.push_back(report::num(local[k], "%.3g"));
            tab.line(row);
        }
        // Harness sanity: the correctly rounded result is within half an ulp of y.
        checks.expect(local[nf - 1] <= 0.5, "CR within 0.5 ulp (" + name + ", " + fp::type_name<T>() + ")");
    }
}

int main() {
    std::printf("t04: pointwise error by position t in the interval\n");
    report::Checks checks;
    run<float>(checks);
    run<double>(checks);
    return checks.exit_code();
}
