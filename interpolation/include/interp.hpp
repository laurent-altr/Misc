// The linear interpolation formulas under test. All of them compute the
// value at x of the segment (x0,y0)-(x1,y1), assuming x0 < x1 and x in [x0,x1].
#pragma once
#include <array>
#include <cmath>

namespace interp {

// A: anchored on the left bound.
template <class T> T left(T x0, T x1, T y0, T y1, T x) {
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

// A': left bound, slope computed first (as when slopes are precomputed).
template <class T> T left_slope(T x0, T x1, T y0, T y1, T x) {
    return y0 + (x - x0) * ((y1 - y0) / (x1 - x0));
}

// A'': left bound, relative position t computed first.
template <class T> T left_t(T x0, T x1, T y0, T y1, T x) {
    return y0 + ((x - x0) / (x1 - x0)) * (y1 - y0);
}

// B: anchored on the right bound.
template <class T> T right(T x0, T x1, T y0, T y1, T x) {
    return y1 - (x1 - x) * (y1 - y0) / (x1 - x0);
}

// C: barycentric weights.
template <class T> T weights(T x0, T x1, T y0, T y1, T x) {
    T t = (x - x0) / (x1 - x0);
    return (1 - t) * y0 + t * y1;
}

// D: fused multiply-add on the last step.
template <class T> T fma_t(T x0, T x1, T y0, T y1, T x) {
    T t = (x - x0) / (x1 - x0);
    return std::fma(t, y1 - y0, y0);
}

// E: anchored on the nearest bound.
template <class T> T nearest(T x0, T x1, T y0, T y1, T x) {
    return (x - x0 < x1 - x) ? left(x0, x1, y0, y1, x) : right(x0, x1, y0, y1, x);
}

// F: C++20 std::lerp (exact at the ends, monotonic).
template <class T> T std_lerp(T x0, T x1, T y0, T y1, T x) {
    return std::lerp(y0, y1, (x - x0) / (x1 - x0));
}

template <class T> using Fn = T (*)(T, T, T, T, T);

template <class T> struct Formula {
    const char* id;
    const char* name;
    Fn<T> f;
};

template <class T> std::array<Formula<T>, 8> formulas() {
    return {{
        {"A",   "left",       left<T>},
        {"A'",  "left_slope", left_slope<T>},
        {"A''", "left_t",     left_t<T>},
        {"B",   "right",      right<T>},
        {"C",   "weights",    weights<T>},
        {"D",   "fma",        fma_t<T>},
        {"E",   "nearest",    nearest<T>},
        {"F",   "std::lerp",  std_lerp<T>},
    }};
}

}  // namespace interp
