# Dynamic FORMAT Plan Experiment Log

This log records the major approaches tried while developing the dynamic FORMAT
parser, the result of each approach, and what survived into the current design.

## Starting Point

The initial problem was that exact, hand-written FORMAT kernels were much faster
than the dynamic implementation, but exact kernels were too brittle.  They only
matched a few complete FORMAT strings, such as:

- `GT:AB:AD:DP:GQ:PL`
- `GT:AD:DP:GQ:PL`
- `GT:AB:AD:DP:GQ:PGT:PID:PL`
- `GT:AD:DP:GQ:PGT:PID:PL`

The target became: recognize useful structure at the FORMAT-tag level, remain
general across subsets/supersets/reordered tags, and fall back to production
htslib whenever the optimized parser could not prove byte-identical output.

## Exact CCDG Kernels

The first high-performance path was a set of exact kernels for dominant CCDG
likelihood layouts.  They proved the upper-bound target: on the 10k CCDG subset,
exact mode was roughly 1.6 s user versus 2.6 s baseline.

Result: useful as a performance oracle, but removed from the production
candidate because exact string matching did not satisfy the generality goal.

## Dynamic Likelihood Shape Executor

Next, the parser used header/type/order information to recognize a likelihood
shape rather than exact tag names:

```text
GT2, optional FLOAT1, INT[n_allele], INT1, INT1,
optional STR1, optional STR1, INT[n_allele * (n_allele + 1) / 2]
```

This was selected by type/order/width rather than names such as `AD` and `PL`.
It validated allele count, observed vector counts, GT syntax, separators, sample
count, and phase-string widths per row.

Result: it closed much of the performance gap.  On one 10k CCDG run, dynamic
shape was within about 6% of exact user time while remaining byte-identical.

Why it did not survive: it reintroduced a shape-specific executor family.  That
was useful evidence, but the MVP goal shifted toward one composable per-tag
executor before adding any generation/specialization layer.

## Cached Shape Classification

The dynamic shape attempt initially paid repeated failed probes on non-likelihood
workloads.  Caching deterministic shape facts per `(header, FORMAT)` plan fixed
that.  The full 1000G GT-only workload stopped paying over a million failed
likelihood-shape probes.

Result: retained as a lesson for future specialization.  The current composable
plan still caches by `(header, FORMAT)`.

## GT-Only Fast Path

A tiny `FORMAT=GT` / diploid `GT2` executor was added and gave a large speedup
on the full 1000G chr22 genotype VCF, cutting dynamic-mode user time from about
9.1 s to about 5.6 s in that intermediate architecture.

Result: the direct `GT2` insight survived, but not as a separate GT-only
executor.  The current composable executor direct-writes leading `GT2` rows when
safe.

## Integer Parse And Encode Tightening

Several low-risk parser/encoder refinements were tried:

- fixed-width integer vector parsers for common AD/PL widths;
- positive integer fast path before falling back to full signed/missing parsing;
- integer range tracking with a `has_special` bit so int8/int16 encoding can skip
  sentinel checks only when the parser proved no missing/vector-end values.

Result: retained.  These fit the generic per-op architecture and helped recover
some likelihood-heavy performance.

## Likelihood Row-Op Elision

In the shape-executor phase, row-op construction was removed from the dynamic
likelihood strict path so the executor could consume cached plan indices and
row-local widths directly.

Result: useful for the old shape executor, but not retained once the MVP pivoted
to the composable row-op model.

## Composable MVP Pivot

The architecture pivoted to:

```text
FORMAT/header -> per-tag compiled ops -> one composable executor -> fallback
```

The dynamic path stopped routing through separate GT-only, likelihood-shape,
fixed-numeric, and measured-general executor ladders.  Instead, it builds one
row-local op list from header metadata and parses supported ops in FORMAT order.

Result: retained.  This is the current design because it supports tag-level
composition for rows such as `GT:AD`, `GT:AD:DP:XX:PL`, reordered fields, and
supersets with normal header-described tags.

Tradeoff: broader composability lost some of the microkernel speed from the
likelihood shape executor.

## Production Hardening

Several hardening passes made the composable MVP safer and faster:

- tightened `GT` compile validation to require `Type=String,Number=1`;
- added malformed-but-readable `GT` header coverage;
- restored direct writes for leading fixed-encoding ops (`GT2`, `FLOAT1`);
- routed generic `INTN` widths 4, 6, and 10 through fixed-width counted parsers;
- removed unused dynamic likelihood-shape scaffolding;
- added underfilled vector compaction for fixed-width vector fields.

Result: retained.  The dynamic path became broader and reduced unnecessary
whole-row fallback while preserving byte-identical output.

## Reverted Or Removed Work

Removed:

- exact CCDG kernels;
- dynamic likelihood-shape executor scaffolding;
- optional SIMD tab-scanning front-end;
- shape-stat benchmark plumbing;
- legacy `exact`/`interp` timing rows in the benchmark harness.

Tested and reverted:

- pointer-increment / reduced-bookkeeping hot-loop rewrite.  It stayed
  byte-correct but slowed targeted likelihood-heavy benchmarks.

## Dynamic-Only Production Trim

After removing exact and SIMD paths, the optimized entry became:

```text
HTS_VCF_FORMAT_PLAN enabled -> dynamic per-tag plan -> composable executor -> production fallback
```

`HTS_VCF_FORMAT_PLAN=1`, `interp`, and `general` now route through the same
dynamic executor.  Benchmarks label only `HTS_VCF_FORMAT_PLAN=1` as `plan`.

Large-corpus post-trim user-time highlights:

| Input | Baseline | Plan | Result |
|---|---:|---:|---|
| CCDG 10k | 2.62 s | 2.25 s | faster, partial fallback |
| 1000G chr22 full GT | 26.05 s | 7.98 s | major win |
| Large CCDG-like synthetic | 4.24 s | 3.78 s | modest win |
| Large float/string | 2.93 s | 2.97 s | near parity/slightly slower |
| Two-string float negative | 2.28 s | 2.56 s | slower |

## bcftools Production Check

A clean bcftools `develop` worktree was built against this htslib branch and run
with `bcftools view --no-version -Ob -l 0 [--threads N]`.

All planned outputs compared byte-identical to baseline.

| Input | Threads | Baseline real | Plan real | Speedup |
|---|---:|---:|---:|---:|
| 1000G chr22 full GT | 0 | 27.48 s | 8.99 s | 3.06x |
| 1000G chr22 full GT | 4 | 26.71 s | 6.94 s | 3.85x |
| Large CCDG-like synthetic | 0 | 4.43 s | 3.94 s | 1.12x |
| Large CCDG-like synthetic | 4 | 3.47 s | 3.02 s | 1.15x |

## Main Lessons

- Tag-level composition is the right MVP boundary; exact full FORMAT strings are
  too brittle.
- Whole-row fallback keeps correctness manageable, but makes one unsupported tag
  enough to lose the optimized path.
- Sample-rich GT-only VCFs are the clearest production win.
- Likelihood-heavy workloads benefit, but generic per-op dispatch and string /
  measured-width handling still leave performance on the table.
- Future executor generation or shape-specialized families may be worth adding
  after the composable MVP is stable.
