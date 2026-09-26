// Reference values in __float128 (113-bit significand). The inputs are the
// stored T values, so the reference is the exact segment between the stored
// points, up to ~2^-113 relative error.
#pragma once
#include <quadmath.h>
#include <string>

#include "fp_utils.hpp"

namespace ref {

using Q = __float128;

template <class T> Q value(T x0, T x1, T y0, T y1, T x) {
    Q qx0 = x0, qx1 = x1, qy0 = y0, qy1 = y1, qx = x;
    Q t = (qx - qx0) / (qx1 - qx0);
    // Blend from the nearest end so the reference is accurate everywhere.
    return t <= Q(0.5) ? qy0 + t * (qy1 - qy0) : qy1 - (1 - t) * (qy1 - qy0);
}

// Correctly rounded result in T: the best any formula could return.
template <class T> T rounded(T x0, T x1, T y0, T y1, T x) {
    return static_cast<T>(value(x0, x1, y0, y1, x));
}

// Error of `got` against the exact value, in ulps of the exact value
// (measured in T).
template <class T> double err_ulps(T got, Q exact) {
    T e = static_cast<T>(exact);
    T u = fp::ulp(e);
    return static_cast<double>((Q(got) - exact) / Q(u));
}

}  // namespace ref
