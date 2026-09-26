# Linear interpolation — floating-point sensitivity experiments

Goal: measure how the **floating-point result** of a linear interpolation
`y(x)` on a table `(X, Y)` behaves when `x` moves by one ulp (ε), and how that
depends on **how the formula is written** (anchored on the left bound, the
right bound, weights, fma, ...).
Mathematical interpolation error is out of scope: the "truth" is the exact
linear segment between the stored points `(x0,y0)`, `(x1,y1)`.

## 1. Formulas to compare

For `x` in `[x0, x1]`, `dx = x1 - x0`, `dy = y1 - y0`:

| id  | name              | expression                                  | notes |
|-----|-------------------|---------------------------------------------|-------|
| A   | left anchor       | `y0 + (x - x0) * dy / dx`                   | exact at `x0` |
| A'  | left, slope first | `y0 + (x - x0) * (dy / dx)`                 | slope may be precomputed per interval |
| A'' | left, t first     | `y0 + ((x - x0) / dx) * dy`                 | |
| B   | right anchor      | `y1 - (x1 - x) * dy / dx`                   | exact at `x1` |
| C   | weights           | `t = (x-x0)/dx; (1-t)*y0 + t*y1`            | may fail `y0 == y1 ⇒ y == y0` |
| D   | fma               | `std::fma(t, dy, y0)`                       | single rounding on the last step |
| E   | nearest anchor    | A if `t < 0.5` else B                       | best accuracy, possible jump at the switch |
| F   | `std::lerp`       | `std::lerp(y0, y1, t)` (C++20)              | libstdc++ guarantees exactness at ends, monotonicity |

Each is a `template <class T>` function so `float`, `double`, `long double`
can all be tested.

## 2. Properties measured

1. **Endpoint exactness**: `f(x0) == y0`, `f(x1) == y1` (bitwise).
2. **Constant preservation**: `y0 == y1 ⇒ f(x) == y0` for all `x`.
3. **Monotonicity under ε**: if `dy > 0`, `f(nextafter(x, +inf)) >= f(x)`.
   Count violations (and their size in ulps of `y`).
4. **ε-discrepancy**: `Δf = f(x⁺) - f(x)` compared with the exact
   `Δy_ref = (x⁺ - x) * dy/dx`. Report in ulps of `y`:
   max, mean, histogram. This is the core question.
5. **Accuracy vs position**: ulp error of `f(x)` against a reference,
   bucketed by `t ∈ [0,1]` (10 bins). Expectation: A is best near `t=0`,
   B near `t=1` (the product `(x-xa)*slope` has an absolute error proportional
   to its magnitude, i.e. to the distance to the anchor).
6. **Continuity at knots**: at an interior knot `x_k`, evaluate with interval
   `[k-1,k]` and `[k,k+1]`, and at `nextbefore/nextafter(x_k)`; report jumps.
7. **Anchor switch (formula E)**: jump / non-monotonicity around `t = 0.5`.

Reference value: compute in `__float128` (gcc, `-lquadmath`) from the same
stored `x0,x1,y0,y1,x`, then round to `T`. For `float` tests, `double` is
already enough but `__float128` keeps one code path.

## 3. Data sets (`cases.hpp`)

Chosen to expose different floating-point regimes:

| case | x0, x1            | y0, y1              | why |
|------|-------------------|---------------------|-----|
| unit | 1, 2              | 0, 1                | baseline, `x - x0` exact (Sterbenz) |
| off  | 1e6, 1e6+1        | 3, 7                | large `x`, few bits for position (timestamps) |
| yoff | 0, 1              | 1e8, 1e8+1          | large `y` offset, tiny slope relative to `y` |
| sign | 0, 1              | -1, 1               | result crosses 0: relative error blows up near 0 |
| flat | 0, 1              | 0.1, 0.1            | constant preservation |
| steep| 1, 1+2^-20        | 0, 1e6              | tiny interval, huge slope: ε in x = big step in y |
| neg  | -3, 5             | 2, -7               | decreasing, mixed signs |
| rand | random tables     | random              | statistical summary |

## 4. Sampling strategy

- `float` on small intervals: **exhaustive** walk `x → nextafter(x)` over
  every representable `x` in `[x0, x1]` (e.g. `[1,2]` = 2^23 values, < 1 s).
- `double`: exhaustive walks over windows of ~10^6 ulps placed at `t≈0`,
  `t≈0.5`, `t≈1`, plus ~10^7 uniform random `x` (fixed seed for
  reproducibility).

## 5. Directory layout

```
interpolation/
  PLAN.md
  Makefile
  include/
    interp.hpp      // formulas A..F, templated on T
    fp_utils.hpp    // ulp(), ulp_distance() via bit_cast, next_up/down, to_string of hex floats
    reference.hpp   // __float128 reference
    cases.hpp       // data sets of section 3
    report.hpp      // tiny table/CSV printer, CHECK macro (no external framework)
  tests/
    t01_endpoints.cpp      // properties 1, 2
    t02_monotonic.cpp      // property 3
    t03_eps_discrepancy.cpp// property 4 (core)
    t04_accuracy_by_t.cpp  // property 5
    t05_knots.cpp          // property 6
    t06_anchor_switch.cpp  // property 7
  results/                 // CSV output (git-ignored or committed snapshots)
  plot.py                  // optional: plots from CSV
```

Each test is a standalone `main()`, prints a human-readable table per
(case × formula × type) and writes a CSV. Tests that check *guaranteed*
properties (A exact at x0, B at x1, F monotonic) return non-zero on failure;
the others are purely observational.

## 6. Build matrix (Makefile)

Floating-point results depend on compiler flags, so the same sources are
built in several configurations into `build/<cfg>/`:

| cfg      | flags |
|----------|-------|
| strict   | `-std=c++20 -O2 -ffp-contract=off` |
| contract | `-std=c++20 -O2 -march=native -ffp-contract=fast` (gcc may emit FMA silently) |
| fast     | `-std=c++20 -O2 -ffast-math` |
| O0       | `-std=c++20 -O0 -ffp-contract=off` |
| x87      | `-std=c++20 -O2 -m32 -mfpmath=387` (excess precision, if multilib available) |

`make` builds & runs `strict`; `make matrix` runs all and diffs the CSVs.
Always pass `-ffp-contract` explicitly: gcc's default differs between
`-std=c++20` and `-std=gnu++20`.

## 7. Expected outcomes (hypotheses to confirm or refute)

- A exact at `x0`, not always at `x1`; B the mirror image.
- ε-discrepancy of A grows with `t`, of B with `1-t`; D (fma) and E reduce it.
- C breaks constant preservation and can be non-monotonic.
- E has the best pointwise accuracy but can show a non-monotonic step at `t=0.5`.
- F (`std::lerp`) monotonic and exact at ends, at some cost in speed.
- `off` case: the position `x - x0` is quantized coarsely; `Δf` is either 0
  or a large jump, independent of formula (input quantization dominates).
- `-ffp-contract=fast` / `-ffast-math` change results of A/B/C bitwise.

## 8. Steps

1. Scaffolding: directory, Makefile, `fp_utils.hpp`, `report.hpp`.
2. `interp.hpp` + `reference.hpp` + `t01_endpoints` (sanity of the harness).
3. `t02_monotonic`, `t03_eps_discrepancy` — the core experiment.
4. `t04_accuracy_by_t`, `t05_knots`, `t06_anchor_switch`.
5. Build matrix, CSV output, optional `plot.py`.
6. Write `README.md` summarizing observed results vs section 7.
