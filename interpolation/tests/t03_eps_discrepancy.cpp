// Property 4 (core experiment): how the computed step f(x+) - f(x), with
// x+ = next_up(x), compares with the exact step of the segment, and with the
// step of the correctly rounded result ("CR", the best achievable).
// Steps are measured in ulps of the segment scale max(|y0|,|y1|): near a zero
// crossing the ulp of y itself becomes tiny and would dominate. The t<.1 and
// t>.9 columns show whether the behaviour depends on the distance to the bound
// the formula is anchored on.
#include <cmath>
#include <cstdio>

#include "cases.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "reference.hpp"
#include "report.hpp"
#include "sampling.hpp"

struct Acc {
    long n = 0, ne_cr = 0;
    double sum = 0, max = 0;
    void add(double err, bool differs) {
        ++n;
        ne_cr += differs;
        sum += err;
        max = std::max(max, err);
    }
    double mean() const { return n ? sum / n : 0; }
    double pct() const { return n ? 100.0 * ne_cr / n : 0; }
};

struct Stat { Acc all, lo, hi; };  // all x, t < 0.1, t > 0.9

template <class T> void run() {
    std::printf("\n=============== %s ===============\n", fp::type_name<T>());
    auto segs = cases::fixed<T>();
    auto rnd = cases::random<T>(16);
    for (auto& s : rnd) s.name = "rand";
    segs.insert(segs.end(), rnd.begin(), rnd.end());

    auto fms = interp::formulas<T>();
    const std::size_t nf = fms.size() + 1;  // + CR baseline
    report::Csv csv(std::string("t03_eps_discrepancy_") + fp::type_name<T>(),
                    "case,formula,range,steps,max_err_ulps,mean_err_ulps,pct_differs_from_cr");

    sampling::Budget budget;
    std::size_t i = 0;
    while (i < segs.size()) {
        const std::string name = segs[i].name;
        std::vector<Stat> st(nf);
        for (; i < segs.size() && segs[i].name == name; ++i) {
            const auto& s = segs[i];
            const ref::Q u = fp::ulp(std::max(std::fabs(s.y0), std::fabs(s.y1)));
            sampling::steps(s, budget, [&](T x) {
                T xp = fp::next_up(x);
                ref::Q e0 = ref::value(s.x0, s.x1, s.y0, s.y1, x);
                ref::Q e1 = ref::value(s.x0, s.x1, s.y0, s.y1, xp);
                ref::Q d_exact = e1 - e0;
                T c0 = static_cast<T>(e0), c1 = static_cast<T>(e1);
                ref::Q d_cr = ref::Q(c1) - ref::Q(c0);
                double t = static_cast<double>((ref::Q(x) - s.x0) / (ref::Q(s.x1) - s.x0));
                for (std::size_t k = 0; k < nf; ++k) {
                    ref::Q d = d_cr;
                    if (k < fms.size())
                        d = ref::Q(fms[k].f(s.x0, s.x1, s.y0, s.y1, xp)) - ref::Q(fms[k].f(s.x0, s.x1, s.y0, s.y1, x));
                    double err = std::fabs(static_cast<double>((d - d_exact) / u));
                    bool differs = d != d_cr;
                    st[k].all.add(err, differs);
                    if (t < 0.1) st[k].lo.add(err, differs);
                    if (t > 0.9) st[k].hi.add(err, differs);
                }
            });
        }

        std::printf("\n-- case %s  (%ld steps)\n", name.c_str(), st[0].all.n);
        report::Table tab({"formula", "max err", "mean err", "%!=CR", "max t<.1", "max t>.9", "%!=CR t<.1", "%!=CR t>.9"}, 10);
        for (std::size_t k = 0; k < nf; ++k) {
            const std::string id = k < fms.size() ? fms[k].id : "CR";
            const std::string fname = k < fms.size() ? fms[k].name : "correctly_rounded";
            const Stat& s = st[k];
            tab.line({id, report::num(s.all.max), report::num(s.all.mean()), report::num(s.all.pct(), "%.2f"),
                      s.lo.n ? report::num(s.lo.max) : "-", s.hi.n ? report::num(s.hi.max) : "-",
                      s.lo.n ? report::num(s.lo.pct(), "%.2f") : "-", s.hi.n ? report::num(s.hi.pct(), "%.2f") : "-"});
            for (auto [range, a] : {std::pair{"all", &s.all}, {"t<0.1", &s.lo}, {"t>0.9", &s.hi}})
                csv.row(name, fname, range, a->n, a->max, a->mean(), a->pct());
        }
    }
}

int main() {
    std::printf("t03: discrepancy of the step f(next_up(x)) - f(x)\n");
    std::printf("err     = |computed step - exact step| in ulps of max(|y0|,|y1|) (CR row: best achievable)\n");
    std::printf("%%!=CR   = %% of steps differing from the step of the correctly rounded result\n");
    run<float>();
    run<double>();
    return 0;
}
