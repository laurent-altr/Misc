# Linear interpolation: floating-point sensitivity experiments

This directory compares ways of writing a linear interpolation `y(x)` on a
segment `(x0,y0)-(x1,y1)`. The comparison is purely about floating-point
behaviour: the "truth" is the exact segment through the stored points, and
approximating the underlying function is out of scope. The main question is
how the result moves when `x` moves by one ulp (ε), and whether that depends
on the bound the formula is anchored on. See [PLAN.md](PLAN.md) for the
design.

```
make              # build + run the strict configuration (-O2 -ffp-contract=off)
make run CFG=fast # one configuration: strict, contract, fast, O0 (x87 if -m32 works)
make matrix       # every configuration, then list the CSVs that differ from strict
make clean
```

This requires g++ with C++20 and libquadmath. The reference values are
computed in `__float128`. Each test prints tables to stdout. The output is
also saved as `results/<cfg>/<test>.txt` and as CSV files. A test exits
non-zero only when a *guaranteed* property is violated. `make matrix` takes
about 7 minutes.

## Formulas

Every formula below is written for `x` in `[x0,x1]`, with `dx = x1-x0` and
`dy = y1-y0`.

| id  | formula |
|-----|---------|
| A   | `y0 + (x-x0)*dy/dx` (left bound) |
| A'  | `y0 + (x-x0)*(dy/dx)` (precomputed slope) |
| A'' | `y0 + ((x-x0)/dx)*dy` (t first) |
| B   | `y1 - (x1-x)*dy/dx` (right bound) |
| C   | `(1-t)*y0 + t*y1` |
| D   | `fma(t, dy, y0)` |
| E   | A when `x-x0 < x1-x`, else B (nearest bound) |
| F   | `std::lerp(y0, y1, t)` |
| CR  | correctly rounded exact value (the best achievable; baseline only) |

## Tests

| test | what it measures |
|------|------------------|
| `t01_endpoints` | `f(x0)==y0`, `f(x1)==y1`, constants preserved, result within `[y0,y1]` |
| `t02_monotonic` | `f(next_up(x))` never goes backwards |
| `t03_eps_discrepancy` | **core**: the step `f(next_up(x)) - f(x)` compared with the exact step, split by `t<0.1` / `t>0.9` |
| `t04_accuracy_by_t` | pointwise error in 10 bins of `t` |
| `t05_knots` | continuity and monotonicity at the interior knots of random tables |
| `t06_anchor_switch` | monotonicity around E's switch point at mid-interval |

Errors are expressed in ulps of `max(|y0|,|y1|)`. Near a zero crossing, the
ulp of `y` itself becomes tiny, so relative errors explode for every formula
that is not correctly rounded. The `ulp(y)` column of t04 shows this.

## Results (gcc 13.3, x86-64)

### 1. The answer depends on which bound the formula is anchored on

From t03, `rand` case, double, strict build. The table shows the percentage of
one-ulp steps whose size differs from the step of the correctly rounded result:

| formula | near x0 (t<0.1) | near x1 (t>0.9) |
|---------|-----------------|-----------------|
| A (left)  | **0.97 %** | 57.5 % |
| B (right) | 44.8 % | **0.92 %** |
| E (nearest) | 0.97 % | 0.92 % |
| C (weights) | 44.1 % | 49.6 % |
| F (lerp) | 18.0 % | 48.3 % |

Near its anchor, a formula computes `y_anchor + small`, and the small
correction carries only a tiny absolute error. Far from its anchor, the
product `(x - x_anchor)*slope` is large, and its rounding error (about one
ulp of the product) lands on the result. Both the pointwise error (t04: A goes
from 0.76 to 2.7 ulps across the interval, and B does the reverse) and the
ε-step error follow the distance to the anchor. The `neg` and `odd` cases
show the same mirror pattern, up to 100 % of steps being off near the "wrong"
bound.

### 2. Consequences at the bounds and knots

- A is always exact at `x0` but not at `x1`, and B is the mirror image. In
  t01 (double), A misses `y1` on 16 % of random segments, and in float on
  42 %.
- In a table (t05), the value just before a knot therefore comes from a
  formula that may miss `Y[k]`. A causes a discontinuity of up to 3 ulps at
  about 8 % of the knots (double) and a few monotonicity breaks. Only E, F
  and C are continuous at the knots.

### 3. Monotonicity

- Every formula built only from correctly rounded operations that are each
  monotone in `x` is monotone: A, A', A'', B, D and F. No violation was found
  in any segment, and t02 checks this.
- C (weights) is not monotone. It is the sum of a decreasing term and an
  increasing term, and it also breaks `y0==y1 ⇒ y==y0`: up to 12 % of the
  samples on flat segments return a different value.
- E (nearest bound) is the most accurate formula overall, but it can step
  backwards at the switch point. t06 found this on 13 of 20 000 float
  segments (by up to 12 ulps), and on 4 double segments once FMA
  contraction is enabled.
- F (`std::lerp`) is monotone, exact at both ends, and continuous at the
  knots, but in the middle of the interval its accuracy is closer to C's.

### 4. Compiler flags change the results

- `-O0` gives bit-identical results to `-O2 -ffp-contract=off`, as expected on
  x86-64 with SSE2 arithmetic. The x87 build was not available on the test
  machine (no 32-bit multilib).
- `-march=native -ffp-contract=fast` lets gcc fuse `a*b+c` into an FMA
  silently. This changed every test:
  - A is exact at `x1` on only 647 of 1008 segments instead of 846.
  - Knot discontinuities double (A: 27 637 knots instead of 11 981).
  - E becomes non-monotone in double.
  - A' becomes much better, since it effectively turns into
    `fma(x-x0, slope, y0)`: 17 % of its steps differ from CR instead of 28 %,
    and its max error drops to 1 ulp.
- `-ffast-math` changed a few double results slightly (reassociation). No
  guaranteed property broke in these runs, but none is guaranteed under it.

### 5. Regimes with no difference

`unit` ([1,2]→[0,1]), `off` (large `x`, 128 representable points in the
interval) and `steep` (1024 points) are exact for every formula: the operands
are powers of two or exactly representable, so every intermediate result is
exact. With a large `y` offset (`yoff`), every formula except C returns the
correctly rounded value. There the quantization of `y` dominates everything
else.

## Practical takeaways

- To get exact values at the bounds of an interval, use E (nearest bound) or
  F (`std::lerp`).
- For monotonicity, avoid E and C. A, B, D and F are safe under strict
  IEEE evaluation.
- For both properties, use F. Patching A with `if (x == x1) return y1` is not
  enough: A can overshoot `y1` just before `x1` (t01, "out range"), so the
  patch can break monotonicity.
- Always set `-ffp-contract` explicitly when reproducibility matters.
