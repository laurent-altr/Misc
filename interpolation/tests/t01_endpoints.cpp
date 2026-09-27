// Properties 1 and 2: exactness at the bounds, constant preservation, and
// results staying inside [min(y0,y1), max(y0,y1)].
#include <algorithm>
#include <cstdio>

#include "cases.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "report.hpp"
#include "sampling.hpp"

template <class T> void run(report::Checks& checks) {
    std::printf("\n=== %s ===\n", fp::type_name<T>());
    auto segs = cases::fixed<T>();
    auto rnd = cases::random<T>(1000);
    segs.insert(segs.end(), rnd.begin(), rnd.end());

    // Flat segments: y0 == y1 with various values.
    std::vector<cases::Segment<T>> flats;
    for (const auto& s : rnd) flats.push_back({s.name + "_flat", s.x0, s.x1, s.y0, s.y0});

    report::Csv csv(std::string("t01_endpoints_") + fp::type_name<T>(),
                    "formula,exact_x0,exact_x1,segments,flat_ok,flat_samples,out_of_range,range_samples");
    report::Table tab({"formula", "exact@x0", "exact@x1", "flat ok", "out range"}, 13);

    for (const auto& fm : interp::formulas<T>()) {
        long ex0 = 0, ex1 = 0, flat_ok = 0, flat_n = 0, out = 0, range_n = 0;
        for (const auto& s : segs) {
            ex0 += fm.f(s.x0, s.x1, s.y0, s.y1, s.x0) == s.y0;
            ex1 += fm.f(s.x0, s.x1, s.y0, s.y1, s.x1) == s.y1;
            T lo = std::min(s.y0, s.y1), hi = std::max(s.y0, s.y1);
            sampling::uniform(s, 2000, [&](T x) {
                T y = fm.f(s.x0, s.x1, s.y0, s.y1, x);
                out += !(y >= lo && y <= hi);
                ++range_n;
            });
        }
        for (const auto& s : flats) {
            sampling::uniform(s, 2000, [&](T x) {
                flat_ok += fm.f(s.x0, s.x1, s.y0, s.y1, x) == s.y0;
                ++flat_n;
            });
        }
        long n = static_cast<long>(segs.size());
        tab.line({fm.id, std::to_string(ex0) + "/" + std::to_string(n), std::to_string(ex1) + "/" + std::to_string(n),
                  report::num(100.0 * flat_ok / flat_n, "%.2f%%"), std::to_string(out)});
        csv.row(fm.name, ex0, ex1, n, flat_ok, flat_n, out, range_n);

        const std::string id = fm.id;
        const std::string who = id + " (" + fp::type_name<T>() + ")";
        if (id == "A" || id == "A'" || id == "A''" || id == "E" || id == "F" || id == "G" || id == "H" || id == "H'" || id == "L" || id == "W") checks.expect(ex0 == n, who + " exact at x0");
        if (id == "B" || id == "E" || id == "F" || id == "G" || id == "H" || id == "H'" || id == "W") checks.expect(ex1 == n, who + " exact at x1");
        if (id != "C") checks.expect(flat_ok == flat_n, who + " preserves constants");
    }
}

int main() {
    std::printf("t01: exactness at bounds, constant preservation, range\n");
    std::printf("exact@x0/x1: segments where f(x0)==y0 / f(x1)==y1 bitwise\n");
    std::printf("flat ok: fraction of samples on y0==y1 segments returning exactly y0\n");
    std::printf("out range: samples outside [min(y0,y1), max(y0,y1)]\n");
    report::Checks checks;
    run<float>(checks);
    run<double>(checks);
    return checks.exit_code();
}
