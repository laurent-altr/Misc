// Small floating-point helpers: ulp, neighbours, ordered integer mapping.
#pragma once
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <type_traits>

namespace fp {

template <class T> T next_up(T x) { return std::nextafter(x, std::numeric_limits<T>::infinity()); }
template <class T> T next_down(T x) { return std::nextafter(x, -std::numeric_limits<T>::infinity()); }

// Size of the gap above |x| (the ulp of x); for 0 the smallest subnormal.
template <class T> T ulp(T x) {
    x = std::fabs(x);
    return next_up(x) - x;
}

template <class T> struct bits;
template <> struct bits<float>  { using type = std::int32_t; };
template <> struct bits<double> { using type = std::int64_t; };

// Maps floats to integers so that consecutive floats map to consecutive
// integers (+0 and -0 both map to 0).
template <class T> std::int64_t ordered(T x) {
    using I = typename bits<T>::type;
    I i = std::bit_cast<I>(x);
    return i < 0 ? -static_cast<std::int64_t>(i & std::numeric_limits<I>::max()) : i;
}

// Number of representable values between a and b (signed, b - a),
// saturated to the int64 range (double intervals crossing 0 can overflow).
template <class T> std::int64_t ulp_distance(T a, T b) {
    __int128 d = static_cast<__int128>(ordered(b)) - ordered(a);
    constexpr std::int64_t m = std::numeric_limits<std::int64_t>::max();
    return d > m ? m : d < -m ? -m : static_cast<std::int64_t>(d);
}

template <class T> const char* type_name() {
    if constexpr (std::is_same_v<T, float>) return "float";
    else return "double";
}

template <class T> std::string hex(T x) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%a", static_cast<double>(x));
    return buf;
}

template <class T> std::string dec(T x) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", std::numeric_limits<T>::max_digits10, static_cast<double>(x));
    return buf;
}

}  // namespace fp
