// The linear interpolation formulas under test. All of them compute the
// value at x of the segment (x0,y0)-(x1,y1), assuming x0 < x1 and x in [x0,x1].
#pragma once
#include <array>
#include <cmath>
#include <type_traits>

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

// G: nearest bound, one fused multiply-add: fma(x - xa, slope, ya).
// The term added to ya is at most |dy|/2 and is rounded only once.
template <class T> T nearest_fma(T x0, T x1, T y0, T y1, T x) {
    const T s = (y1 - y0) / (x1 - x0);
    return (x - x0 < x1 - x) ? std::fma(x - x0, s, y0) : std::fma(x - x1, s, y1);
}

// Error-free transformations need every operation rounded separately. With
// FMA contraction enabled (-ffp-contract=fast and an FMA-capable -march), gcc
// fuses a product into the following addition, so the error term computed from
// the rounded product no longer matches. This attribute disables contraction
// for the function and everything inlined into it. It does NOT protect against
// -ffast-math: that needs a separate translation unit compiled without it.
#if defined(__GNUC__) && !defined(__clang__)
#define INTERP_NO_CONTRACT __attribute__((optimize("fp-contract=off")))
#else
#define INTERP_NO_CONTRACT
#endif

// Error-free transformation: a + b == s + e exactly.
template <class T> void two_sum(T a, T b, T& s, T& e) {
    s = a + b;
    T bb = s - a;
    e = (a - (s - bb)) + (b - bb);
}

// H: nearest bound, compensated. Keeps the rounding errors of dx, dy, of the
// slope, of x - xa (bits of x lost when |x - xa| >> |x|) and of the product,
// and adds them back before the final rounding (double-length accuracy).
template <class T> INTERP_NO_CONTRACT T nearest_comp(T x0, T x1, T y0, T y1, T x) {
    const bool left_side = x - x0 < x1 - x;
    const T xa = left_side ? x0 : x1, ya = left_side ? y0 : y1;
    T dx, dx_e, dy, dy_e;
    two_sum(x1, -x0, dx, dx_e);
    two_sum(y1, -y0, dy, dy_e);
    const T s = dy / dx;
    // (dy + dy_e) / (dx + dx_e) ~= s + s_lo
    const T s_lo = (std::fma(-s, dx, dy) + dy_e - s * dx_e) / dx;
    T h, h_e;
    two_sum(x, -xa, h, h_e);                            // x - xa == h + h_e exactly
    const T p = h * s;
    const T p_e = std::fma(h, s, -p);                   // h * s == p + p_e exactly
    T sum, sum_e;
    two_sum(ya, p, sum, sum_e);
    return sum + (sum_e + (p_e + (h * s_lo + h_e * s)));
}

// Error-free product without fma (Dekker): a * b == p + e exactly, using
// Veltkamp's split of each factor into two half-width parts. Valid while
// |a|, |b| stay far from the overflow threshold.
template <class T> void split(T a, T& hi, T& lo) {
    constexpr T factor = std::is_same_v<T, float> ? T(4097) : T(134217729);  // 2^ceil(p/2) + 1
    const T c = factor * a;
    hi = c - (c - a);
    lo = a - hi;
}

template <class T> void two_prod(T a, T b, T& p, T& e) {
    p = a * b;
    T ah, al, bh, bl;
    split(a, ah, al);
    split(b, bh, bl);
    e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
}

// H': same as H without fma: the exact products use two_prod instead.
template <class T> INTERP_NO_CONTRACT T nearest_comp_nofma(T x0, T x1, T y0, T y1, T x) {
    const bool left_side = x - x0 < x1 - x;
    const T xa = left_side ? x0 : x1, ya = left_side ? y0 : y1;
    T dx, dx_e, dy, dy_e;
    two_sum(x1, -x0, dx, dx_e);
    two_sum(y1, -y0, dy, dy_e);
    const T s = dy / dx;
    T q, q_e;
    two_prod(s, dx, q, q_e);                            // s * dx == q + q_e
    const T r = (dy - q) - q_e;                         // dy - s * dx, exact
    const T s_lo = (r + dy_e - s * dx_e) / dx;
    T h, h_e;
    two_sum(x, -xa, h, h_e);
    T p, p_e;
    two_prod(h, s, p, p_e);
    T sum, sum_e;
    two_sum(ya, p, sum, sum_e);
    return sum + (sum_e + (p_e + (h * s_lo + h_e * s)));
}

template <class T> using wide_t = std::conditional_t<std::is_same_v<T, float>, double, long double>;

// W: nearest bound computed in a wider type (double for float, x87 80-bit
// long double for double), rounded once at the end.
template <class T> T nearest_wide(T x0, T x1, T y0, T y1, T x) {
    using W = wide_t<T>;
    return static_cast<T>(nearest<W>(W(x0), W(x1), W(y0), W(y1), W(x)));
}

template <class T> using Fn = T (*)(T, T, T, T, T);

template <class T> struct Formula {
    const char* id;
    const char* name;
    Fn<T> f;
};

template <class T> std::array<Formula<T>, 12> formulas() {
    return {{
        {"A",   "left",       left<T>},
        {"A'",  "left_slope", left_slope<T>},
        {"A''", "left_t",     left_t<T>},
        {"B",   "right",      right<T>},
        {"C",   "weights",    weights<T>},
        {"D",   "fma",        fma_t<T>},
        {"E",   "nearest",    nearest<T>},
        {"F",   "std::lerp",  std_lerp<T>},
        {"G",   "nearest_fma",  nearest_fma<T>},
        {"H",   "nearest_comp", nearest_comp<T>},
        {"H'",  "nearest_comp_nofma", nearest_comp_nofma<T>},
        {"W",   "nearest_wide", nearest_wide<T>},
    }};
}

}  // namespace interp
