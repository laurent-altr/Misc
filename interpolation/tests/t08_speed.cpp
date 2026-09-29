// Rough cost of each formula: ns per call over random segments and x,
// called through a function pointer (same overhead for every formula).
// Observational only; numbers depend on the machine and the flags.
#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

#include "fp_utils.hpp"
#include "interp.hpp"
#include "report.hpp"

template <class T> void run() {
    std::printf("\n=== %s ===\n", fp::type_name<T>());
    constexpr int N = 1 << 20, R = 20;
    std::mt19937_64 g(1);
    std::uniform_real_distribution<double> u(0, 1), v(-100, 100);
    struct In { T x0, x1, y0, y1, x; };
    std::vector<In> in(N);
    for (auto& e : in) {
        e.x0 = T(v(g));
        e.x1 = T(e.x0 + 0.001 + 10 * u(g));
        e.y0 = T(v(g));
        e.y1 = T(v(g));
        e.x = T(e.x0 + u(g) * (e.x1 - e.x0));
    }
    report::Table tab({"formula", "ns/call"}, 10);
    report::Csv csv(std::string("t08_speed_") + fp::type_name<T>(), "formula,ns_per_call");
    for (const auto& fm : interp::formulas<T>()) {
        volatile T sink = 0;
        double best = 1e30;
        for (int r = 0; r < R; ++r) {
            T acc = 0;
            auto t0 = std::chrono::steady_clock::now();
            for (const auto& e : in) acc += fm.f(e.x0, e.x1, e.y0, e.y1, e.x);
            auto t1 = std::chrono::steady_clock::now();
            sink = sink + acc;
            best = std::min(best, std::chrono::duration<double, std::nano>(t1 - t0).count() / N);
        }
        tab.line({fm.id, report::num(best, "%.2f")});
        csv.row(fm.name, best);
    }
}

int main() {
    std::printf("t08: cost per call (best of 20 runs over 2^20 calls)\n");
    run<float>();
    run<double>();
    return 0;
}
