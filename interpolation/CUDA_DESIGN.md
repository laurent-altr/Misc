# CUDA version: design notes

This document describes a design only; there is no code yet. It covers a CUDA
port of the compensated interpolation L (see [README.md](README.md)) for the
target application: a finite element solver that interpolates strain/stress
curves, with **one interpolation per thread** (one integration point per
thread). The value `x` of each point changes slowly between solver
iterations.

The CPU conclusions do not all carry over. On the CPU, the best layout was a
per-lane cache of the interval data. On a GPU it is most likely the opposite:
keep one small shared table resident in the cache and store only the interval
index per point. The reasons follow.

---

## 1. What each thread needs

Per call, thread `i` handles one integration point:

| data | size (double) | access pattern |
|------|---------------|----------------|
| `x[i]` (strain) | 8 B | read, coalesced |
| `k[i]` (current interval, the initial guess) | 4 B | read, coalesced; written only when it changes |
| `y[i]` (stress) | 8 B | written, coalesced |
| interval data for `k[i]`: `x0, y0, s, s_lo` | 32 B | **gather**: depends on `k[i]` |
| upper bound `x1 = X[k+1]` (interval check) | 8 B | gather, next to the previous one |

On the GPU, `fma` is native and exact, so the per-interval data of L needs
only 4 values:
- `x0`, the anchor;
- `y0`;
- `s`, the rounded slope;
- `s_lo`, the correction of the slope.

The Veltkamp split used by the fma-free H' is not needed. The per-call
arithmetic is about 15 FP64 operations: two TwoSums, one product with its
error computed by `fma`, and the final sum.

The kernel is therefore **memory-bound**. The design is about bytes moved,
not about arithmetic.

---

## 2. Where to put the table (Y values and interval data)

### Recommended: global memory, read-only path, record of 32 bytes per interval

- Store the table as an **array of 32-byte records**, one per interval:
  `{x0, y0, s, s_lo}` in double (16 bytes `{x0, y0, s, s_lo}` in float),
  aligned on 32 bytes.
  - One record is exactly one 32-byte memory sector, so each thread's gather
    costs a single sector transaction, loaded as two 128-bit loads
    (`double2`), or one `float4` in float.
  - A structure of arrays (separate `X[]`, `Y[]`, `S[]`, `S_lo[]`) would
    touch 4 sectors per thread for the same data. That layout was best on the
    CPU (vector gathers per array). On the GPU, the array of records is
    better.
  - The upper bound `x1` of interval `k` is the `x0` of record `k+1`, which
    is the adjacent sector and usually in the same 128-byte line. Keep one
    extra sentinel record at the end, whose `x0` is the table's last knot.
- Read it through the **read-only data path**: `const __restrict__` pointer
  parameters, or `__ldg`, so the loads can use the L1/texture cache.
- A curve of 1 000 knots is 32 KB in double. It stays in **L2** (several MB
  on current GPUs) across kernel launches and solver iterations, and largely
  in L1 during a kernel. The gathers are then served from cache, not DRAM.
- Upload the tables once, keep them resident on the device for the whole
  simulation, and re-upload only when a curve changes.

### Alternatives, and why not by default

| option | verdict |
|--------|---------|
| **`__constant__` memory** (64 KB) | The constant cache broadcasts only when all threads of a warp read the **same address**. Different intervals within a warp are serialized, up to 32 times. It is only good if all points of a warp are almost always in the same interval. That can happen when neighbouring elements have similar strain, but it is not a safe default. |
| **Shared memory** (copy the table per block at kernel start) | Every block must copy the whole table (32 KB per block for 1 000 double knots) to do one interpolation per thread. The copy costs more than the work unless the table is very small (a few hundred bytes) or each thread does many interpolations. Random 8-byte accesses also cause bank conflicts. Reconsider it only if profiling shows L1/L2 misses on the table. |
| **Texture objects with hardware linear filtering** | **Do not use for this purpose.** The hardware computes the interpolation weight with about 8 fractional bits (fixed point), which is far worse than even the naive formula and produces exactly the staircase this work tries to avoid. Textures also do not filter doubles. Reading through a texture without filtering brings nothing over the read-only path on current GPUs. |
| **Per-thread cache of the interval data** (the CPU winner) | Each call would read 40-48 bytes of per-point state from DRAM (coalesced, but not cacheable: millions of points exceed L2), instead of 4 bytes of `k` plus a gather that hits L2. On the GPU this moves about 3-4 times more DRAM bytes. Keep only `k` per point. |

### Several curves (several materials)

- Concatenate all tables in one device buffer, and give each material an
  offset (`first record`, `number of records`). A point then needs its
  material's offset: store it per element, or better, in a small table per
  material indexed by the element's material id.
- **Order the points so that a warp handles one material**, which is the
  usual element ordering by material. The warp's gathers then stay inside one
  small table, which is good for L1 hits. Mixing materials within a warp is
  correct but less cache-friendly.
- If all tables together exceed a few MB (hundreds of curves with thousands
  of knots), measure L2 hit rates before deciding anything else.

---

## 3. Per-point state and the interval search

- Keep `k[i]` (int32) in global memory, in point order: coalesced read, and
  written back **only if it changed**, which saves a store in almost every
  call. This is the GPU form of the "initial guess". No flag for the side
  (left or right bound) is needed, because L always anchors on `x0`.
- Search: starting from `k`, step down while `x < x0[k]`, step up while
  `x >= x0[k+1]`, clamped to the table.
  - With slowly changing `x`, the loop body almost never runs, so warp
    divergence is negligible.
  - Bound the walk to a few steps and fall back to a binary search. This
    protects the first iteration, where the guesses are poor, and any large
    jump in strain.
- First call: initialize `k` with a binary search, either in a separate
  kernel or with the same kernel using a bad guess.
- Points outside the table extrapolate with the first or last interval, the
  same convention as the CPU code. Use finite sentinels, not infinities, if
  fast-math flags are ever used.

---

## 4. Thread and memory mapping

- One thread per integration point, one-dimensional grid, and a grid-stride
  loop so the kernel works for any number of points. Blocks of 128-256
  threads.
- `x`, `k` and `y` in structure-of-arrays layout indexed by point
  (`x[i]`, not `element[i].x`), so that a warp reads 32 consecutive values:
  256 bytes, fully coalesced.
- If the solver stores strains inside a per-element structure (array of
  structures), the interpolation kernel will read with a stride. Prefer
  separate contiguous arrays for the values this kernel reads and writes, or
  fuse the interpolation into the kernel that already has `x` in registers.
- **Fusing** into the solver's existing constitutive kernel is probably the
  biggest real-world gain: `x` is already in a register, and `y` is consumed
  immediately. The interpolation then costs only the `k` read and the table
  gather.

---

## 5. Numerical correctness on the GPU

The CPU lessons apply, with CUDA-specific switches:

- **FMA contraction is on by default in nvcc** (`--fmad=true`). The compiler
  may fuse a multiply and an add, which breaks the error-free
  transformations exactly as seen on the CPU (t09). Two options:
  - write the error-free parts with the intrinsics that are **never
    contracted**: `__dadd_rn`, `__dsub_rn`, `__dmul_rn` (and `__fadd_rn`,
    `__fsub_rn`, `__fmul_rn` in float), plus an explicit `fma()` for the
    product error. This is the recommended option: correct whatever the
    compilation flags;
  - or compile that translation unit with `--fmad=false`. This is simpler,
    but it slows down all other code in that unit.
- **`--use_fast_math`**: it affects single precision (approximate division,
  flush-to-zero). Avoid it for this code in float. At least keep
  `-ftz=false`, because flushing tiny intermediate values to zero breaks
  TwoSum near 0, and keep `-prec-div=true`. Double precision operations stay
  IEEE.
- **Precompute the table on the host** with the existing CPU code
  (`make_comp_seg` / `chunk::make_table`: slope and slope correction), then
  copy it. The device then does no division at all, and the table is
  bit-identical to the CPU one. CUDA's `+`, `-`, `*` and `fma` are IEEE
  round-to-nearest, so **the GPU results must be bit-identical to the CPU
  `left_comp_nofma`**. This is the main acceptance test.

---

## 6. Double or float: check the GPU model

- FP64 throughput depends strongly on the GPU:
  - data-center GPUs (A100, H100 and similar) run FP64 at half the FP32 rate;
  - consumer and workstation GPUs (GeForce, most RTX) run it 32-64× slower.
- On a data-center GPU the kernel stays memory-bound in double. On a
  consumer GPU, the roughly 15 FP64 operations may become the bottleneck.
  Measure before optimizing anything else.
- In float, L is correctly rounded in float (99.9 % of steps equal to the
  correctly rounded step, t07). If the solver can run the interpolation in
  float, a record is 16 bytes (`float4`, one load) and all the arithmetic
  runs at full speed.

---

## 7. Expected cost (to verify)

Per point and per call, in double:
- DRAM traffic: 8 B read (`x`) + 4 B read (`k`) + 8 B written (`y`), plus
  rarely 4 B written (`k`). That is about **20 bytes**.
- Cache traffic: one or two 32-byte sectors for the table, from L1/L2.

At 1-3 TB/s of DRAM bandwidth, that is on the order of 0.01-0.02 ns per point
when the kernel is memory-bound. This is two orders of magnitude below the
CPU cost (about 0.7-1 ns). Kernel launch overhead (a few µs) dominates unless
there are at least about 10^5-10^6 points per launch, which is another
argument for fusing the interpolation into an existing kernel.

The naive formula would move the same bytes. On a GPU with good FP64
throughput, **L should cost about the same as the naive formula**, since the
extra arithmetic is hidden behind the memory traffic.

---

## 8. Validation and benchmark plan

1. **Bit-exactness:** random tables and random `x`, including the special
   probes of t07 (near `x≈0`, near `y≈0`, near the knots, outside the
   table). The GPU `y` must equal the CPU `left_comp_nofma` bit for bit, in
   float and in double.
2. **Interval tracking:** after each call, `X[k] <= x < X[k+1]` for every
   point (with the clamping at the ends), under slow walks, fast walks and
   jumps.
3. **Flag robustness:** build with and without `--fmad=false` and
   `--use_fast_math`. Results must not change when the error-free parts use
   the `_rn` intrinsics.
4. **Performance** (Nsight Compute), measuring:
   - effective bandwidth against the device peak;
   - L1/L2 hit rate of the table gathers;
   - warp divergence of the search loop;
   - the cost relative to the naive formula, in the same kernel.
5. Variants to compare:
   - table in global memory (read-only path) against a shared-memory copy;
   - 32-byte records against a structure of arrays;
   - float against double;
   - standalone kernel against fused into the constitutive kernel.

---

## 9. Open questions (they change the recommendations)

- **GPU model:** is FP64 fast (data-center) or slow (consumer)?
- **Number and size of the curves:** knots per curve, and number of
  materials. This decides whether all tables fit in L2 comfortably.
- **Points per call:** below about 10^5 per launch, overhead dominates and
  fusion becomes essential.
- **Precision:** is float acceptable for the interpolated stress, or is
  double required?
- **Data layout in the solver:** are strains and stresses already in
  contiguous per-point arrays, or inside per-element structures?
