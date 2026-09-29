// Main goal: when x runs over consecutive representable values, y should
// follow the correctly rounded result (CR) as closely as possible, with no
// extra "stairs" (flat runs followed by jumps). CR itself has stairs when a
// one-ulp move of x changes y by less than one ulp of y; those cannot be
// avoided. Everything is measured in ulps of y at x (local ulps).
//
// Probes: runs of N consecutive x starting near t=0, t=0.5, t=1, near x=0
// (intervals containing 0), near y=0 (segments crossing 0), and at random t.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <random>

#include "cases.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "reference.hpp"
#include "report.hpp"

template <class T> struct Probe {
    cases::Segment<T> s;
    T start;
    std::string where;
};

template <class T> std::vector<Probe<T>> probes() {
    auto segs = cases::fixed<T>();
    // y = 2x through the origin: the interval contains x = 0 and y = 0.
    segs.push_back({"cross", T(-1), T(3), T(-2), T(6)});
    std::vector<Probe<T>> out;
    auto add_fixed = [&](const cases::Segment<T>& s) {
        const T dx = s.x1 - s.x0;
        for (auto [t, w] : {std::pair{0.02, "t~0"}, {0.5, "t~.5"}, {0.98, "t~1"}})
            out.push_back({s, T(s.x0 + T(t) * dx), w});
        if (s.x0 < 0 && s.x1 > 0) out.push_back({s, std::ldexp(T(1), -10), "x~0"});
        if ((s.y0 < 0) != (s.y1 < 0)) {
            T xz = T(s.x0 - s.y0 * (dx / (s.y1 - s.y0)));
            for (int i = 0; i < 64; ++i) xz = fp::next_down(xz);
            out.push_back({s, xz, "y~0"});
        }
    };
    for (const auto& s : segs) add_fixed(s);
    std::mt19937_64 g(99);
    std::uniform_real_distribution<double> u(0, 1);
    for (const auto& s : cases::random<T>(300, 555)) {
        out.push_back({s, T(s.x0 + u(g) * (s.x1 - s.x0)), "rand"});
        add_fixed(s);
    }
    return out;
}

struct Stat {
    long steps = 0, same = 0, stalls = 0, backwards = 0;
    std::int64_t max_dev = 0, extra_flat = 0;
};

template <class T> void run() {
    std::printf("\n=============== %s ===============\n", fp::type_name<T>());
    constexpr int N = 2048;
    auto fms = interp::formulas<T>();
    const std::size_t nf = fms.size() + 1;
    const std::vector<std::string> wheres{"t~0", "t~.5", "t~1", "x~0", "y~0", "rand"};
    std::map<std::string, std::vector<Stat>> st;
    for (const auto& w : wheres) st[w].assign(nf, {});
    std::vector<Stat> total(nf);

    // Consecutive y values of the demo probes, printed as steps in ulps.
    struct Demo { std::string title; std::vector<std::vector<std::int64_t>> d; };
    std::vector<Demo> demos;

    for (const auto& pr : probes<T>()) {
        const auto& s = pr.s;
        std::vector<T> xs;
        for (T x = pr.start; x < s.x1 && xs.size() < N + 1; x = fp::next_up(x))
            if (x >= s.x0) xs.push_back(x);
        if (xs.size() < 2) continue;
        std::vector<T> cr(xs.size());
        for (std::size_t i = 0; i < xs.size(); ++i) cr[i] = ref::rounded(s.x0, s.x1, s.y0, s.y1, xs[i]);
        const bool inc = s.y1 >= s.y0;
        auto longest_flat = [](const std::vector<T>& y) {
            std::int64_t run = 0, best = 0;
            for (std::size_t i = 1; i < y.size(); ++i) best = std::max(best, run = y[i] == y[i - 1] ? run + 1 : 0);
            return best;
        };
        const std::int64_t cr_flat = longest_flat(cr);

        const bool demo = (s.name == "neg" && pr.where == "t~1") || (s.name == "cross" && pr.where == "x~0") ||
                          (s.name == "odd" && pr.where == "t~.5") || (s.name == "yoff" && pr.where == "t~.5");
        Demo dm{s.name + " " + pr.where + " (x from " + fp::dec(xs[0]) + ")", {}};

        for (std::size_t k = 0; k < nf; ++k) {
            std::vector<T> y(xs.size());
            for (std::size_t i = 0; i < xs.size(); ++i)
                y[i] = k < fms.size() ? fms[k].f(s.x0, s.x1, s.y0, s.y1, xs[i]) : cr[i];
            Stat a;
            a.extra_flat = longest_flat(y) - cr_flat;
            std::vector<std::int64_t> d;
            for (std::size_t i = 0; i < xs.size(); ++i) {
                a.max_dev = std::max(a.max_dev, std::abs(fp::ulp_distance(cr[i], y[i])));
                if (i == 0) continue;
                std::int64_t df = fp::ulp_distance(y[i - 1], y[i]), dc = fp::ulp_distance(cr[i - 1], cr[i]);
                ++a.steps;
                a.same += df == dc;
                a.stalls += df == 0 && dc != 0;
                a.backwards += inc ? df < 0 : df > 0;
                if (demo && i <= 32) d.push_back(df);
            }
            for (Stat* t : {&st[pr.where][k], &total[k]}) {
                t->steps += a.steps;
                t->same += a.same;
                t->stalls += a.stalls;
                t->backwards += a.backwards;
                t->max_dev = std::max(t->max_dev, a.max_dev);
                t->extra_flat = std::max(t->extra_flat, a.extra_flat);
            }
            if (demo) dm.d.push_back(d);
        }
        if (demo) demos.push_back(dm);
    }

    auto id = [&](std::size_t k) { return std::string(k < fms.size() ? fms[k].id : "CR"); };
    auto name = [&](std::size_t k) { return std::string(k < fms.size() ? fms[k].name : "correctly_rounded"); };

    std::printf("\nSteps equal to the CR step (%%), by probe position:\n");
    std::vector<std::string> head{"formula"};
    head.insert(head.end(), wheres.begin(), wheres.end());
    head.push_back("all");
    report::Table t1(head, 8);
    report::Csv csv(std::string("t07_stairs_") + fp::type_name<T>(),
                    "formula,where,steps,same_as_cr,stalls,backwards,max_dev_ulps,extra_flat");
    for (std::size_t k = 0; k < nf; ++k) {
        std::vector<std::string> row{id(k)};
        for (const auto& w : wheres) {
            const Stat& a = st[w][k];
            row.push_back(report::num(100.0 * a.same / a.steps, "%.2f"));
            csv.row(name(k), w, a.steps, a.same, a.stalls, a.backwards, a.max_dev, a.extra_flat);
        }
        const Stat& a = total[k];
        row.push_back(report::num(100.0 * a.same / a.steps, "%.2f"));
        csv.row(name(k), "all", a.steps, a.same, a.stalls, a.backwards, a.max_dev, a.extra_flat);
        t1.line(row);
    }

    std::printf("\nAll probes (%ld steps): stalls = flat steps where CR moves; dev = |y - CR| in ulps;\n"
                "extra flat = longest flat run minus CR's longest flat run on the same probe\n", total[0].steps);
    report::Table t2({"formula", "%=CR", "stalls", "backwards", "max dev", "extra flat"}, 12);
    for (std::size_t k = 0; k < nf; ++k) {
        const Stat& a = total[k];
        t2.line({id(k), report::num(100.0 * a.same / a.steps, "%.3f"), std::to_string(a.stalls), std::to_string(a.backwards),
                 report::num(double(a.max_dev)), std::to_string(a.extra_flat)});
    }

    for (const auto& dm : demos) {
        std::printf("\nsuccessive steps of y in ulps, %s:\n", dm.title.c_str());
        for (std::size_t k = 0; k < nf; ++k) {
            std::printf("  %-4s", id(k).c_str());
            for (auto v : dm.d[k]) std::printf("%3lld", static_cast<long long>(v));
            std::printf("\n");
        }
    }
}

int main() {
    std::printf("t07: staircase effect over consecutive representable x\n");
    run<float>();
    run<double>();
    return 0;
}
