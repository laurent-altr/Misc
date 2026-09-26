// Data sets: segments chosen to expose different floating-point regimes.
// Magnitudes depending on the precision are expressed with D = number of
// significand bits of T, so float and double see the same regime.
#pragma once
#include <cmath>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace cases {

template <class T> struct Segment {
    std::string name;
    T x0, x1, y0, y1;
};

template <class T> std::vector<Segment<T>> fixed() {
    constexpr int D = std::numeric_limits<T>::digits;
    // 2^(D-8) has an ulp of 2^-7: an interval of length 1 holds 128 values.
    const T big = std::ldexp(T(1), D - 8);
    // 1 + 1024 ulps of 1.
    const T tiny_dx = std::ldexp(T(1), -(D - 1 - 10));
    return {
        {"unit",  T(1),   T(2),            T(0),   T(1)},
        {"off",   big,    big + 1,         T(3),   T(7)},
        {"yoff",  T(0),   T(1),            big,    big + 1},
        {"sign",  T(0),   T(1),            T(-1),  T(1)},
        {"flat",  T(0),   T(1),            T(0.1), T(0.1)},
        {"steep", T(1),   T(1) + tiny_dx,  T(0),   T(1e6)},
        {"neg",   T(-3),  T(5),            T(2),   T(-7)},
        {"odd",   T(0.1), T(0.7),          T(0.3), T(1.9)},
    };
}

// Random segments, reproducible (fixed seed).
template <class T> std::vector<Segment<T>> random(int n, unsigned seed = 12345) {
    std::mt19937_64 g(seed);
    std::uniform_real_distribution<double> pos(-100, 100), len(-3, 1), val(-100, 100);
    std::vector<Segment<T>> out;
    for (int i = 0; i < n; ++i) {
        T x0 = T(pos(g));
        T x1 = T(x0 + std::pow(10.0, len(g)));
        if (!(x1 > x0)) continue;
        out.push_back({"rand" + std::to_string(i), x0, x1, T(val(g)), T(val(g))});
    }
    return out;
}

// Random sorted tables for the knot tests.
template <class T> struct Table {
    std::vector<T> x, y;
};

template <class T> std::vector<Table<T>> random_tables(int n, int knots, bool monotone,
                                                        unsigned seed = 777) {
    std::mt19937_64 g(seed);
    std::uniform_real_distribution<double> step(1e-3, 3), val(-50, 50), inc(0, 5);
    std::vector<Table<T>> out;
    for (int i = 0; i < n; ++i) {
        Table<T> t;
        double x = val(g), y = val(g);
        for (int k = 0; k < knots; ++k) {
            T xs = T(x);
            if (!t.x.empty() && !(xs > t.x.back())) { x += step(g); continue; }
            t.x.push_back(xs);
            t.y.push_back(T(monotone ? y : val(g)));
            x += step(g);
            y += inc(g);
        }
        out.push_back(std::move(t));
    }
    return out;
}

}  // namespace cases
