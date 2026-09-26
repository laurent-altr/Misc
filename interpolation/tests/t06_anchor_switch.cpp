// Property 7: formula E (and G, H, W) uses the left bound in the first half of the interval
// and the right one in the second half. Around the switch the two formulas
// round differently, so the result may step backwards. This scans a window
// of consecutive x around the midpoint of many random increasing segments.
#include <cstdio>

#include "cases.hpp"
#include "fp_utils.hpp"
#include "interp.hpp"
#include "report.hpp"

template <class T> void run() {
    std::printf("\n=== %s ===\n", fp::type_name<T>());
    auto segs = cases::random<T>(20000, 2024);
    for (auto& s : segs)
        if (s.y1 < s.y0) std::swap(s.y0, s.y1);

    const std::int64_t half = 256;
    report::Table tab({"formula", "segments", "seg w/ viol", "violations", "worst ulps"}, 12);
    report::Csv csv(std::string("t06_anchor_switch_") + fp::type_name<T>(), "formula,segments,segments_with_violation,violations,worst_ulps");
    bool shown = false;
    for (const auto& fm : interp::formulas<T>()) {
        const std::string id = fm.id;
        if (id != "A" && id != "B" && id != "E" && id != "G" && id != "H" && id != "W") continue;
        long bad_segs = 0, viol = 0;
        std::int64_t worst = 0;
        for (const auto& s : segs) {
            T x = s.x0 + (s.x1 - s.x0) / 2;
            for (std::int64_t i = 0; i < half && x > s.x0; ++i) x = fp::next_down(x);
            bool bad = false;
            T prev = fm.f(s.x0, s.x1, s.y0, s.y1, x);
            for (std::int64_t i = 0; i < 2 * half && x < s.x1; ++i) {
                T xn = fp::next_up(x);
                T y = fm.f(s.x0, s.x1, s.y0, s.y1, xn);
                std::int64_t back = fp::ulp_distance(y, prev);
                if (back > 0) {
                    ++viol;
                    bad = true;
                    worst = std::max(worst, back);
                    if (!shown && id == "E") {
                        shown = true;
                        std::printf("example: x0=%s x1=%s y0=%s y1=%s\n", fp::dec(s.x0).c_str(), fp::dec(s.x1).c_str(),
                                    fp::dec(s.y0).c_str(), fp::dec(s.y1).c_str());
                        std::printf("  f(%s) = %s\n  f(%s) = %s  (smaller)\n", fp::hex(x).c_str(), fp::dec(prev).c_str(),
                                    fp::hex(xn).c_str(), fp::dec(y).c_str());
                    }
                }
                prev = y;
                x = xn;
            }
            bad_segs += bad;
        }
        tab.line({id, std::to_string(segs.size()), std::to_string(bad_segs), std::to_string(viol), std::to_string(worst)});
        csv.row(fm.name, segs.size(), bad_segs, viol, worst);
    }
}

int main() {
    std::printf("t06: monotonicity around the mid-interval switch of formula E (vs A and B)\n");
    run<float>();
    run<double>();
    return 0;
}
