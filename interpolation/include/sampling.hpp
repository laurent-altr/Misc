// Choice of the x values at which a segment is probed.
#pragma once
#include <algorithm>
#include <cstdint>
#include <random>

#include "cases.hpp"
#include "fp_utils.hpp"

namespace sampling {

struct Budget {
    std::int64_t exhaustive_max = std::int64_t(1) << 23;  // walk every x if the interval is this small
    std::int64_t window = std::int64_t(1) << 16;          // else: consecutive x near t=0, 0.5, 1
    std::int64_t random = std::int64_t(1) << 16;          // ... plus uniformly random x
};

// Calls f(x) for x in [x0, x1) such that next_up(x) <= x1: every x if the
// interval is small, otherwise three windows of consecutive values plus
// random values. Consecutive calls are consecutive floats inside a window.
template <class T, class F> void steps(const cases::Segment<T>& s, const Budget& b, F&& f) {
    const std::int64_t n = fp::ulp_distance(s.x0, s.x1);
    if (n <= b.exhaustive_max) {
        for (T x = s.x0; x < s.x1; x = fp::next_up(x)) f(x);
        return;
    }
    auto walk = [&](T x, std::int64_t count) {
        for (std::int64_t i = 0; i < count && x < s.x1; ++i, x = fp::next_up(x)) f(x);
    };
    walk(s.x0, b.window);
    T mid = s.x0 + (s.x1 - s.x0) / 2;
    for (std::int64_t i = 0; i < b.window / 2; ++i) mid = fp::next_down(mid);
    walk(mid, b.window);
    T end = s.x1;
    for (std::int64_t i = 0; i < b.window; ++i) end = fp::next_down(end);
    walk(end, b.window);

    std::mt19937_64 g(42);
    std::uniform_real_distribution<double> u(0, 1);
    for (std::int64_t i = 0; i < b.random; ++i) {
        T x = std::clamp(T(s.x0 + u(g) * (s.x1 - s.x0)), s.x0, fp::next_down(s.x1));
        f(x);
    }
}

// Uniformly distributed x in [x0, x1] (every x if the interval is small).
template <class T, class F> void uniform(const cases::Segment<T>& s, std::int64_t count, F&& f) {
    if (fp::ulp_distance(s.x0, s.x1) <= count) {
        for (T x = s.x0; x <= s.x1; x = fp::next_up(x)) f(x);
        return;
    }
    std::mt19937_64 g(43);
    std::uniform_real_distribution<double> u(0, 1);
    for (std::int64_t i = 0; i < count; ++i)
        f(std::clamp(T(s.x0 + u(g) * (s.x1 - s.x0)), s.x0, s.x1));
}

}  // namespace sampling
