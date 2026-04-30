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
HTS_VCF_FORMAT_PLAN enabled -> dynamic per-tag plan -> composable executor -> generic fallback
```

`HTS_VCF_FORMAT_PLAN=1` now routes through the dynamic executor.  Older
`interp` and `general` aliases were later removed during production tightening
so unknown values do not accidentally enable the fast path.

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

## String/Float Shape Boundary

The expanded threaded benchmark exposed two regressions:

- `GT:GL:FT:DP:GQ`
- `GT:FT:PID:GL:DP`

Both schemas were syntactically supported and had zero row-local fallback, but
they were dominated by measured strings plus `Number=G` float vectors.  The
dynamic path had to measure string widths over every sample before parsing, then
still use the general float conversion path, while there were no integer vectors
to amortize that setup.

Result: retained as a conservative support boundary.  The compiler now
negative-caches measured-string plus float-vector schemas that do not also have
integer-vector work, and sends those FORMAT rows to the generic parser.  The
full threaded corpus remained byte-identical.  The two-string float case
improved from a consistent slowdown, roughly 0.86-0.89x, to parity at
1.00-1.01x.  Other integer-heavy likelihood rows stayed on the dynamic path.

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

Result: retained.  The FORMAT-plan tests now compare explicit inclusion and
exclusion sample lists byte-for-byte against production parsing.  A
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

## GIAB and Full CCDG Probe

Four GIAB HG002 VCFs were pulled into `bench/format-shape/large/public/giab`:
NIST v4.2.1 GRCh38 small variants, v5.0q GRCh38 small variants, v5.0q GRCh38
structural variants, and v5.0q CHM13v2.0 small variants.  The bcftools command
suite was run on those files plus the 3,202-sample CCDG 10k slice using
`bench/format-shape/large/bcftools-giab-ccdg-inputs.tsv`.

First result: GIAB v5.0q exposed a real GT correctness bug.  The planned GT2
parser encoded `.|.` as `./.` because missing alleles were stored without the
separator phase bit.  The parser now accepts simple diploid missing/digit
combinations and preserves the phase bit for `.|.`, `0|.`, and `.|0`; an
explicit edge row was added to `test/format-plan-edge.vcf`.

After the fix, all baseline-vs-plan outputs compared `ok`.  Speedups are modest
on GIAB because it is single-sample data: roughly 1.06-1.11x for `view_bcf` and
1.03-1.09x for `query_format`.  CCDG 10k remained in the expected cohort range:
1.13x for `view_bcf`, 1.52x for `query_format`, and 1.10x for `filter_gt`.

The full parent CCDG/1000G high-coverage chr22 VCF was identified as:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz
```

It is 26.0 GiB compressed and requires a local copy for reruns.
The normal command harness materializes complete outputs, which is not practical
for this file: a single `view_bcf -Ob -l 0` baseline output reached 155 GiB
before that run was stopped.  A streaming checksum harness was added so command
outputs can be validated without storing them.

The full CCDG streaming command suite completed with all baseline-vs-plan
checksums comparing `ok`:

| Command | Baseline real | Plan real | Real speedup | Baseline user | Plan user | User speedup |
|---|---:|---:|---:|---:|---:|---:|
| `view_bcf` | 678.46 s | 562.96 s | 1.21x | 476.41 s | 377.47 s | 1.26x |
| `view_sites` | 472.27 s | 403.28 s | 1.17x | 455.70 s | 386.18 s | 1.18x |
| `query_sites` | 71.44 s | 76.78 s | 0.93x | 67.02 s | 72.00 s | 0.93x |
| `query_format` | 124.14 s | 76.88 s | 1.61x | 119.16 s | 72.27 s | 1.65x |
| `stats` | 77.45 s | 77.12 s | 1.00x | 72.86 s | 72.55 s | 1.00x |
| `filter_gt` | 531.20 s | 453.21 s | 1.17x | 512.95 s | 434.35 s | 1.18x |

## Parser Helper Trim

Reviewed the `vcf.c` implementation for duplicated fast-path helper code.  A
first attempt collapsed the fixed-width integer vector parsers into one generic
counted loop.  Correctness held, but likelihood-shaped rows regressed by roughly
10% in the focused benchmark, so that version was rejected.

The retained refactor is intentionally narrower: remove unused non-range integer
vector helpers, remove an unused scalar helper, and centralize only the empty
integer-vector fill case.  The hand-unrolled range parsers for common vector
widths remain because they are part of the measured hot path.

Result: retained.  The final `vcf.c` diff is about 116 fewer deleted helper
lines relative to the previous branch tip, with byte-identical outputs on the
focused GT/likelihood/string corpus.  A repeat likelihood benchmark was neutral:
CCDG-like plan time improved from 4.18 s to 4.12 s, reordered likelihood was
2.66 s to 2.70 s, and multiallelic likelihood improved from 3.01 s to 2.98 s.

## Production Tightening Review

Three review passes focused on production-readiness: code-size risk, correctness
risk, and upstream polish.  The retained implementation changes were deliberately
low-risk:

- `HTS_VCF_FORMAT_PLAN` now enables only on `1`; old `interp` / `general`
  aliases and typo-enables were removed.
- Planner statistics are incremented only when
  `HTS_VCF_FORMAT_PLAN_STATS=1`, avoiding process-global counter writes in
  normal runs.
- The row-op support check was folded into row-op resolution, removing a second
  pass over the FORMAT operation list.
- The row-width bound was made explicit in the planner instead of using an
  inline literal.
- Tests now assert that an unknown value such as `HTS_VCF_FORMAT_PLAN=off`
  behaves like the generic parser.

Result: retained.  `make test/test_view test/test_format_plan_cache`,
the FORMAT-plan parser-output checks, and `test/test_format_plan_cache` pass.
At that point, the `vcf.c` implementation was about 1,594 added lines relative to
`origin/develop`, down from the earlier 1,703-line core.

## Generic Executor Micro-Optimizations

The next pass targeted the generic per-op executor rather than adding new
schema-specific kernels.  Retained changes:

- skip `max_counts` maintenance for row ops that cannot compact;
- update integer min/max directly on the common positive-integer parse path;
- reject over-wide measured `Number=.` / string fields during the measurement
  pass instead of after scanning the full row;
- remove nullable `nread` checks from planner-private integer vector helpers.

Result: retained.  Focused FORMAT tests passed, `git diff --check` was clean,
and the htslib large corpus in
`bench/format-shape/large/results-opt-batch1b` compared byte-identical to
baseline.

| Input | Baseline user | Plan user | User speedup | Hits/fallback |
|---|---:|---:|---:|---:|
| CCDG 10k | 2.50 s | 2.20 s | 1.14x | 8,396 / 1,604 |
| 1000G chr22 full GT | 25.08 s | 8.99 s | 2.79x | 1,103,547 / 0 |
| Large CCDG-like synthetic | 4.02 s | 3.68 s | 1.09x | 20,000 / 0 |
| Large reordered likelihood | 2.91 s | 2.38 s | 1.22x | 20,000 / 0 |
| Large multiallelic likelihood | 3.19 s | 2.64 s | 1.21x | 16,000 / 0 |
| Large float/string | 2.89 s | 2.88 s | 1.00x | 0 / 16,000 |
| Variable phase widths | 2.57 s | 2.44 s | 1.05x | 12,000 / 0 |
| Mixed row-local fallbacks | 2.20 s | 1.83 s | 1.20x | 12,000 / 0 |
| GT-first reordered | 1.73 s | 1.41 s | 1.23x | 12,000 / 0 |
| Two-string float | 2.25 s | 2.24 s | 1.00x | 0 / 12,000 |

One broader structural attempt was rejected: splitting the all-samples loop from
the `keep_samples` loop.  Correctness held, but
`bench/format-shape/large/results-opt-nosubset-split` was slower across the
planned corpus: CCDG 10k plan user time moved from 2.20 s to 2.28 s, 1000G
GT-only from 8.99 s to 9.30 s, and the likelihood-shaped synthetic rows also
regressed.  That change was reverted.

The standard bcftools GIAB/CCDG command corpus was then run against a bcftools
binary explicitly linked to this checkout with:

```sh
make HTSDIR=../htslib-vcf-avx-sanity bcftools
```

All command outputs compared `ok` in
`bench/format-shape/large/results-bcftools-giab-ccdg-opt-batch1`.  The command
profile stayed positive where FORMAT parsing matters and neutral/noisy where it
does not: CCDG 10k `query_format` was 1.55x faster by user time, CCDG 10k
`view_bcf` was 1.12x faster, and GIAB single-sample `query_format` rows were
roughly 1.08-1.12x faster.

## Fallback Reason Counters And Split Width Caps

The next regression investigation focused on CCDG rows that were falling back
because phase-set string fields exceeded the old single planned-width limit.
The implementation now reports fallback reasons under
`HTS_VCF_FORMAT_PLAN_STATS=1`:

- unsupported schema;
- numeric width;
- string width;
- GT shape;
- parse failure;
- separator mismatch;
- sample-count mismatch.

The single width cap was split into a 64-value numeric-vector cap and a
256-byte measured-string cap.  Numeric and string width fallbacks are diagnostic
only: they do not disable a schema that succeeds on nearby rows.

Two string caps were benchmarked.  A 512-byte cap planned all CCDG 10k rows but
had a mixed bcftools-level signal.  The retained 256-byte cap planned 9,861 of
10,000 CCDG rows and left the 139 longest string rows on the generic parser:

```text
vcf-format-plan attempts=10000 hits=9861 fallback=139 parsed_samples=31574922
vcf-format-plan-fallback unsupported=0 numeric_width=0 string_width=139 gt_shape=0 parse=0 separator=0 sample_count=0
```

Result: retained.  Focused tests passed, `git diff --check` was clean, and the
htslib large corpus in `bench/format-shape/large/results-string-cap256-reasons`
compared byte-identical to baseline.  CCDG 10k user time was 2.47 s baseline
versus 2.17 s planned, 1000G chr22 full GT was 24.70 s versus 9.75 s, and the
likelihood-shaped synthetic rows remained faster or neutral.

The bcftools GIAB/CCDG command corpus in
`bench/format-shape/large/results-bcftools-giab-ccdg-cap256` also compared
byte-identical.  CCDG 10k user-time speedups were 1.13x for `view_bcf`, 1.56x
for `query_format`, and 1.10x for `filter_gt`; GIAB FORMAT-query rows were
1.07-1.15x faster, while site-only controls and `stats` remained neutral/noisy.

## Repo Test Harness Hardening

The final hardening pass moved the important small-case checks from the
benchmark directory into the actual htslib test harness.  The bespoke shell test
was removed; the production-facing checks now live in `test/test.pl` as
`test_vcf_format_plan`, while `test/test_format_plan_cache` remains the focused
cache-generation check.

Retained changes:

- successful planned rows call `vcf_parse_format_check7()`, matching the generic
  parser's final FORMAT cardinality validation;
- fallback diagnostics are test-only hooks with `*_for_test` names and are
  emitted only when `HTS_VCF_FORMAT_PLAN_STATS=1`;
- `test_vcf_format_plan` compares planned output against generic output
  byte-for-byte, including selected-sample cases and disabled-control values
  such as `HTS_VCF_FORMAT_PLAN=off`;
- new fixtures cover rollback after partial planned parsing, malformed
  unselected samples under `bcf_hdr_set_samples()`, repeated wide GT values, and
  malformed sample-count failures;
- row-local width fallbacks remain record-local so sparse over-cap string rows
  do not poison CCDG-like schemas.

Result: retained.  `make check` passed with 377/377 tests.  `make
maintainer-check` was attempted but failed before the whitespace/copyright
checks because the local build invoked the C compiler on `test/usepublic.cpp`
with `-std=gnu23`.  The relevant whitespace check and `git diff --check` passed
separately.

The htslib large corpus run written locally under
`bench/format-shape/large/results-prod-hardening2` compared byte-identical to
baseline.  The generated result files are ignored, so the recorded summary is:
CCDG 10k held the expected 9,861 / 139 hit/fallback split, and 1000G chr22 full
GT remained the largest win at 24.61 s baseline user time versus 9.48 s planned.

The latest bcftools GIAB/CCDG command corpus in
`bench/format-shape/large/results-bcftools-giab-ccdg-prod-hardening` also
compared byte-identical.  CCDG 10k user-time speedups were 1.14x for `view_bcf`,
1.56x for `query_format`, and 1.12x for `filter_gt`; GIAB single-sample FORMAT
rows remained modestly positive, as expected.

## Runtime Cooldown Removal

The per-plan runtime cooldown was removed after an A/B pass showed no practical
benefit on realistic workloads.  The cooldown had paused a supported cached
schema after repeated row-local fallbacks, but standard corpus hit/fallback
counts were identical with and without it.  The remaining protection is simpler:
compile-time unsupported schemas are negative-cached, unsupported mixed
string/float shapes are rejected at compile time, and row-local misses fall back
only for that record.

The final no-cooldown parser corpus in
`bench/format-shape/large/results-no-cooldown-final` compared byte-identical to
baseline.  Representative planned user times:

| Input | Baseline user | Planned user | Hits / fallback |
|---|---:|---:|---:|
| CCDG 10k | 2.46 s | 2.16 s | 9,861 / 139 |
| 1000G chr22 full GT | 24.50 s | 9.34 s | 1,103,547 / 0 |
| Large reordered likelihood | 2.89 s | 2.42 s | 20,000 / 0 |
| Large float/string negative | 2.88 s | 2.86 s | 0 / 16,000 |
| Mixed row-local fallbacks | 2.14 s | 1.83 s | 12,000 / 0 |
| Two-string float negative | 2.21 s | 2.22 s | 0 / 12,000 |

Result: retained.  Focused planner tests passed, `test/test_format_plan_cache`
passed, all large parser-corpus outputs compared byte-identical, and
`git diff --check` was clean.

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
