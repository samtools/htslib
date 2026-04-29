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
- replaced the original process-global 16-entry FORMAT plan cache with a
  header-owned, generation-aware, dynamically sized cache that stores both
  supported and unsupported compile results.

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

## Header-Owned Cache Hardening

The static FORMAT plan cache was replaced with private `bcf_hdr_aux_t` state.
The hardened cache:

- grows from 16 to 128 entries;
- stores literal FORMAT strings on the heap, so long schemas are no longer
  rejected by the old fixed key buffer;
- caches unsupported compile results to avoid repeated work;
- clears on `bcf_hdr_sync()` and records a private header generation;
- declines fast planning while `h->dirty` is set.

Result: retained.  `test/format-plan-cache.vcf` now asserts 21/21 planned hits
across more than 16 distinct FORMAT schemas, including one long schema.  The
new `test/test_format_plan_cache` helper verifies that a plan compiled before a
header metadata change is not reused after `bcf_hdr_sync()`.  The large corpus
remained byte-identical after the rewrite, with the same broad performance
profile: 1000G chr22 GT user time at 26.06 s baseline versus 7.96 s planned,
and CCDG 10k at 2.55 s baseline versus 2.24 s planned.

## Profitability Gate For String/Float Shapes

The expanded threaded benchmark exposed two regressions:

- `GT:GL:FT:DP:GQ`
- `GT:FT:PID:GL:DP`

Both schemas were syntactically supported and had zero row-local fallback, but
they were dominated by measured strings plus `Number=G` float vectors.  The
dynamic path had to measure string widths over every sample before parsing, then
still use the general float conversion path, while there were no integer vectors
to amortize that setup.

Result: retained.  The compiler now negative-caches these low-profit schemas and
sends only those FORMAT rows to the production parser.  The full threaded corpus
remained byte-identical.  The two-string float case improved from a consistent
slowdown, roughly 0.86-0.89x, to parity at 1.00-1.01x.  Other integer-heavy
likelihood rows stayed on the dynamic path.

## Selected-Sample Support

The planner originally rejected `h->keep_samples` because sample subsetting
changes the relationship between input sample columns and output BCF sample
slots.  That was conservative but would have made the optimized path invisible
for common `bcftools view -s/-S` style workflows.

The executor now treats the input and output counts separately.  It scans
`h->nsamples_ori` columns when `h->keep_samples` is active, skips unselected
columns with the header bitset, writes retained samples densely, and sets
`v->n_sample` to the retained sample count.  The width-measurement pass follows
the same rule, so measured strings and variable numeric widths are based only on
the samples that will be emitted, matching production htslib's selected-sample
behavior.

Result: retained.  `test/test_format_plan.sh` now compares explicit inclusion
and exclusion sample lists byte-for-byte against production parsing.  A
bcftools run selecting the first two samples from every input completed 40/40
byte-identical comparisons.  The 1000G chr22 GT workload still showed a large
real-time win, from 26.51 s to 9.77 s unthreaded and from 25.99 s to 8.84 s at
4 threads; string/float-heavy negative rows remained near parity.

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

## Broader bcftools Command Check

Added `bench/format-shape/scripts/run_bcftools_command_bench.sh` so the branch
can exercise more than `bcftools view`.  The runner currently covers full BCF
conversion, genotype-dropping conversion, site queries, small FORMAT queries,
`stats`, genotype filters, and an opt-in merge benchmark.  Every command runs
once with `HTS_VCF_FORMAT_PLAN=0` and once with `HTS_VCF_FORMAT_PLAN=1`, then
compares outputs with `cmp`.

Result: retained.  All applicable planned outputs compared byte-identical to
baseline.  FORMAT-heavy commands showed the expected gains: 1000G full GT was
2.79x faster for `view_bcf`, 2.98x faster for `view_sites`, 1.94x faster for
`query_format`, and 1.57x faster for `filter_gt`.  CCDG and reordered
likelihood workloads were smaller but positive.  Site-only queries and `stats`
were mostly neutral, with a few small negative rows that remain useful overhead
watchpoints.

`bcftools merge` was tested through the opt-in `merge_self` command against a
smaller manifest to avoid excessive duplicated-sample output.  All planned merge
outputs compared byte-identical to baseline.  Merge was neutral-to-positive:
small 1000G genotype input improved from 0.14 s to 0.10 s, large CCDG
likelihood improved from 4.50 s to 4.33 s, and large float/string remained
unchanged at 2.69 s.

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
