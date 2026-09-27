# Linear interpolation: floating-point sensitivity experiments

This directory compares ways of writing a linear interpolation `y(x)` on a
segment `(x0,y0)-(x1,y1)`. The comparison is purely about floating-point
behaviour: the "truth" is the exact segment through the stored points, and
approximating the underlying function is out of scope. The main question is
how the result moves when `x` moves by one ulp (ε), and whether that depends
on the bound the formula is anchored on. See [PLAN.md](PLAN.md) for the
design.

**Goal: no stairs.** When `y` is printed for successive representable `x`,
it should follow the correctly rounded values, without flat runs followed by
jumps that the exact segment does not have. The answer is formula **H**
(nearest bound + compensated arithmetic, below). For `float` data, formula
**W** (computed in `double`) is even simpler and just as good.

```
make              # build + run the strict configuration (-O2 -ffp-contract=off)
make run CFG=fast # one configuration: strict, native, contract, fast, O0, noinline (x87 if -m32 works)
make matrix       # every configuration, then list the CSVs that differ from strict
make clean
```

This requires g++ with C++20 and libquadmath. The reference values are
computed in `__float128`. Each test prints tables to stdout. The output is
also saved as `results/<cfg>/<test>.txt` and as CSV files. A test exits
non-zero only when a *guaranteed* property is violated. `make matrix` takes
about 12 minutes. The `fast` configuration is expected to fail t09 (see below).

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
| G   | nearest bound, `fma(x-xa, dy/dx, ya)` |
| H   | nearest bound, compensated: the rounding errors of `dx`, `dy`, the slope, `x-xa` and the product are computed exactly (TwoSum / fma) and added back before the last rounding |
| H'  | H without fma: exact products by Dekker's two_prod (Veltkamp split); bit-identical to H |
| W   | E computed in a wider type (`double` for float data, 80-bit `long double` for double data), rounded once |
| CR  | correctly rounded exact value (the best achievable; baseline only) |

## Tests

| test | what it measures |
|------|------------------|
| `t01_endpoints` | `f(x0)==y0`, `f(x1)==y1`, constants preserved, result within `[y0,y1]` |
| `t02_monotonic` | `f(next_up(x))` never goes backwards |
| `t03_eps_discrepancy` | **core**: the step `f(next_up(x)) - f(x)` compared with the exact step, split by `t<0.1` / `t>0.9` |
| `t04_accuracy_by_t` | pointwise error in 10 bins of `t` |
| `t05_knots` | continuity and monotonicity at the interior knots of random tables |
| `t06_anchor_switch` | monotonicity around the switch point of E, G, H, W at mid-interval |
| `t07_stairs` | **goal**: runs of 2048 consecutive `x` near t≈0, t≈0.5, t≈1, x≈0, y≈0 and random t; steps compared with the CR steps, and the step sequences printed for a few cases |
| `t08_speed` | ns per call |
| `t10_naive_vs_comp` | head-to-head A (naive) vs H': quality summary and speed in three usage patterns |
| `t09_error_free` | `two_sum` / `two_prod` exact under the current flags (vs `__float128`); H vs H' bit for bit |

Errors are expressed in ulps of `max(|y0|,|y1|)`. Near a zero crossing, the
ulp of `y` itself becomes tiny, so relative errors explode for every formula
that is not correctly rounded. The `ulp(y)` column of t04 shows this.

## Results (gcc 13.3, x86-64)

### 0. Staircase: which formula follows the ideal steps (t07)

Percentage of one-ulp steps of `x` where the step of `y` equals the step of
the correctly rounded result. "Stalls" counts steps where `y` stays flat
while CR moves. Strict build.

| formula | float %=CR | double %=CR | double: at y≈0 | double stalls | ns/call (double) |
|---------|-----------:|------------:|---------------:|--------------:|-----------------:|
| A (left)            | 53.0 | 53.8 | 1.2 | 15 537 | 2.5 |
| E (nearest)         | 65.3 | 67.3 | 1.2 |  9 759 | 7.7 |
| F (`std::lerp`)     | 49.0 | 48.2 | 1.2 | 16 106 | 14.5 |
| G (nearest + fma)   | 82.4 | 84.7 | 66.0 | 7 954 | 9.3 (7.9 with hardware FMA) |
| **H (compensated)** | **99.9** | **99.93** | **99.4** | **0** | 14.5 (10.0 with hardware FMA) |
| W (wider type)      | **100.0** | 88.4 | 1.8 | 18 | 10.3 |

The same probe, `y` printed for 32 successive `x` (steps in ulps of `y`,
segment `(-3,2)-(5,-7)`, double, `x` near `x1`):

```
A    0 -2  0 -2 -2  0 -2  0 -2  0 -2  0 -2 -2  0 -2  0 -2 ...   <- stairs of 2 ulps
H   -1 -1 -1 -1 -1 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -2 -1 -1 ...
CR  -1 -1 -1 -1 -1 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -2 -1 -1 ...
```

and for `y = 2x` on `[-1,3]` with `x` starting at 2^-10 (CR moves by 1 ulp at
every step). Every formula except H and W stays flat for about 1000
consecutive `x`, then jumps by about 1000 ulps:

```
A..G  0  0  0  0  0  0  0 ...
H, W  1  1  1  1  1  1  1 ...
```

Stairs have three distinct causes:

1. **The rounding grid of y.** If a one-ulp move of `x` changes the exact `y`
   by less than one ulp of `y`, then even CR is a staircase (case `yoff`).
   No formula can avoid this; only a wider type for `y` can.
2. **The term added to the anchor is coarser than y.** In `ya + p`, `p` is
   rounded to its own grid. When `|p| > |y|` (the anchor is far away, or `y`
   crosses 0 between the anchor and `x`), that grid is coarser than the grid
   of `y`, and `y` moves in steps of 2, 4, ... ulps. Anchoring on the nearest
   bound (E, G) halves `|p|` at most, so it does not remove this cause.
3. **`x - xa` drops bits of `x`.** When `|x - xa|` is much larger than `|x|`
   (an interval containing 0, with `x` near 0), consecutive `x` give the same
   `x - xa`, so `y` stays flat for up to 1000 steps and then jumps.

H removes causes 2 and 3. Its error terms are computed exactly with TwoSum
and fma, and the result is rounded only once, which gives double-length
accuracy (about 106 bits for double). The only steps where it still differs
from CR are at `y≈0`, where heavy cancellation leaves too few correct bits
relative to a tiny `y`. W works for the same reason, but an 80-bit
`long double` only adds 11 bits to a double. That is plenty for `float`
data computed in `double`, but not for double data near `y≈0`.

Guarantees and caveats for H:
- It is exact at both bounds (checked in t01) and preserves constants.
- No backward step was observed (t02, t06), but monotonicity is not proven.
- It needs strict IEEE evaluation. **`-ffast-math` destroys it** (t10 under
  `-ffast-math`: H' drops to 74 % correctly rounded): the compiler
  simplifies the TwoSum error terms to 0, and H drops to 69 % of steps equal
  to CR, which is no better than E. `-ffp-contract=fast` does no harm.

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

### H without fma (H'), and what the compiler may break (t09)

H' replaces each `fma` by Dekker's exact product: each factor is split into
two half-width parts (`c = (2^27+1)*a; hi = c - (c - a); lo = a - hi`, with
2^12+1 for float), and the error of `a*b` is rebuilt from the partial
products. In every correct build, H' returns exactly the same bits as H on all
probes, at the same cost without hardware FMA (about 15 ns).

Can `two_sum` / `two_prod` be inlined? **Yes.** Without `-ffast-math`, gcc
never reassociates floating-point operations, and inlining does not change
the result: the `strict`, `O0` and `noinline` builds give identical results,
and all 2×2^21 error-free transformations are exact. What can break them is a
**flag**, not inlining:

| build | two_sum / two_prod exact | H' steps = CR |
|-------|--------------------------|---------------|
| strict, O0, noinline | yes | 99.9 % |
| `-march=native -ffp-contract=fast`, without protection | yes (in isolation) | **78 %, with 195 backward steps** |
| same, with `INTERP_NO_CONTRACT` (current code) | yes | 99.9 % |
| `-ffast-math` | **no**: error terms simplified to 0 | 68 % |

- **FMA contraction.** Each helper stays exact on its own. But once they are
  inlined together, gcc fuses the product `p = h*s` into the next addition
  (`ya + p`), which then uses the *exact* product while the error term assumes
  the *rounded* `p`. The fix is to disable contraction for the function:
  `__attribute__((optimize("fp-contract=off")))`. Everything inlined into it
  inherits this. The header applies it to H and H' through the
  `INTERP_NO_CONTRACT` macro. Compiling with `-ffp-contract=off` works too.
- **`-ffast-math`.** It allows reassociation, so `(a - (s - bb)) + (b - bb)`
  is simplified to 0. No function attribute fixes this: `optimize("no-fast-math")`
  and the individual `no-associative-math`… options were tried with gcc 13,
  and neither works. The only reliable option is to put the function in its
  own `.cpp` file compiled without `-ffast-math` (verified). It is then called
  normally, but it cannot be inlined into fast-math code.
- **x87** (`-m32 -mfpmath=387`) would also break them, because intermediate
  results are kept with extra precision. This could not be tested here. Use
  SSE2 (`-mfpmath=sse`), which is the default on x86-64.
- Veltkamp's split overflows for `|a|` above about 2^996 (double) or 2^115
  (float).

### Naive (A) vs H': quality and speed (t10)

Quality: 2000 random segments, 512 random `x` each, plus runs of 512
consecutive `x` (at a random `t` and near `x1`), plus 1000 random tables.

| metric | float A | float H' | double A | double H' |
|--------|--------:|---------:|---------:|----------:|
| correctly rounded results | 61.9 % | 100.0 % | 64.7 % | 100.0 % |
| mean error (ulps of y) | 2.9 | 0.00002 | 2.6 | 0 |
| max error (ulps of y, near y≈0) | 87 500 | 1 | 141 000 | 0 |
| steps equal to the CR step | 50.3 % | 100.0 % | 50.9 % | 100.0 % |
| stalls (flat while CR moves) | 6 598 | 0 | 6 154 | 0 |
| backward steps | 0 | 0 | 0 | 0 |
| extra flat run | 68 | 0 | 42 | 0 |
| exact at `x1` | 62.3 % | 100 % | 84.4 % | 100 % |
| knots where the left interval misses `Y[k]` | 37.5 % | 0 % | 15.3 % | 0 % |

Speed in ns per call, on a 2.1 GHz Xeon. Each loop is an out-of-line loop
over plain arrays, so gcc may inline and vectorize the formula. Hp is H' with
the slope and its correction precomputed per segment (`make_comp_seg` /
`eval_comp_seg`, same bits as H').

| pattern | build | A | H' | Hp | H'/A |
|---------|-------|--:|---:|---:|-----:|
| independent (new segment and `x` each call) | `-O2` | 1.1 | 11.3 | - | 10× |
| | `-O3 -march=native` | 1.1 | 10.2 | - | 9× |
| sweep (one segment, array of `x`) | `-O2` | 1.1 | 9.0 | 3.8 | 8× |
| | `-O3 -march=native` | 0.54 | 1.04 | 1.12 | **2×** |
| table (1024 knots, binary search + interpolation) | `-O2` | 56 | 65 | 58 | 1.2× |
| | `-O3 -march=native` | 49 | 58 | 65 | 1.2× |

(double; float is similar. `results/<cfg>/t10_naive_vs_comp.txt` has both.)

- **Independent calls:** H' costs about 10 ns more. Most of it is the
  arithmetic: 2 divisions and about 45 other operations, against 1 division
  and 3 operations for A. The unpredictable choice of the nearest bound adds
  about 4 ns: with `x` always in the left half, H' costs 6.6 ns. A branchless
  selection (bit masks) removes the misprediction but does not make it faster
  overall, so it was not kept.
- **Sweep:** at `-O3 -march=native`, gcc vectorizes both formulas (AVX2,
  4 doubles at a time), and H' is only 2× A. At `-O2`, gcc 13 does not inline
  H' (it is too large for the `-O2` heuristics), so each point pays a call and
  no vectorization. Use `-O3`, or force inlining.
- **Table:** the binary search (about 55 ns with random queries) dominates, so
  H' adds only about 10 ns (+20 %). Precomputing the segment data (Hp) removes
  the divisions but makes no reliable difference here. It pays off only in
  scalar sweeps at `-O2` (3.8 against 9 ns). Timings in this pattern vary by ±10 ns between runs.
- The `INTERP_NO_CONTRACT` attribute does not change speed: with or without
  it, the times are identical. It can be disabled with
  `-DINTERP_NO_CONTRACT=` in a build that already uses `-ffp-contract=off`.

## Practical takeaways

- **For smooth `y` over successive `x` (no extra stairs), use H.** For `float`
  data, W (compute in `double`, round once) is simpler and just as good. Both
  cost about 3-5× a plain formula (roughly 10-15 ns instead of 2.5 ns).
- Cost of H' against the naive formula: about +10 ns per isolated call, 2×
  in a vectorized loop at `-O3 -march=native`, and +20 % when a table search
  is involved.
- Without hardware FMA, use H' (same bits as H). Keep `INTERP_NO_CONTRACT` on
  it so a later `-march=native` build cannot silently break it.
- Never compile H or H' with `-ffast-math`. If needed, isolate it in a translation
  unit compiled with `-fno-fast-math`.
- To get exact values at the bounds of an interval, use E (nearest bound) or
  F (`std::lerp`).
- For monotonicity, avoid E and C. A, B, D and F are safe under strict
  IEEE evaluation.
- For both properties, use F. Patching A with `if (x == x1) return y1` is not
  enough: A can overshoot `y1` just before `x1` (t01, "out range"), so the
  patch can break monotonicity.
- Always set `-ffp-contract` explicitly when reproducibility matters.
