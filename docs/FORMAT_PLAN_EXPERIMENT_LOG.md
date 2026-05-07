# FORMAT Planner Experiment Log

This branch is an experiment in making VCF FORMAT parsing faster while keeping
the default HTSlib behavior unchanged.  The planner remains disabled unless
`HTS_VCF_FORMAT_PLAN=1` is set.

## Current Result

The current implementation uses a dynamic per-tag planner with a conservative
whole-column fallback contract.  It does not rely on width-specific integer
micro-specializations.  The executor now uses one `INTVEC` row kind for integer
vectors instead of separate `INT2`, `INT3`, and selected fixed-width parsers.

The simplification removed about 283 lines from `vcf.c`:

```text
vcf.c | 300 ++++--------------------------------------------------------------
1 file changed, 17 insertions(+), 283 deletions(-)
```

## Correctness Evidence

Commands run after the simplification:

```sh
make test/test_view test/test_format_plan_cache bgzip tabix
./test/test_format_plan_cache
perl test/test.pl -F test_vcf_format_plan
git diff --check
```

The Perl FORMAT-plan fragment passed 21/21 tests.  A streaming `cmp` comparison
between the pre-simplification planner and simplified planner also passed for
all 27 local corpus inputs.

## Local Corpus Simplification Benchmark

Reference implementation: detached worktree at
`0999351c040793ae9763b2f9ba7550b45fd9bd11`.

Compared variants:

- `before`: planner before integer-vector simplification.
- `simplified`: `codex/simplify-format-plan-executor`.

Both variants used:

```sh
HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 test/test_view -b -l 0 INPUT >/dev/null
```

Each input had one warmup per variant followed by three paired measured runs.
The input corpus was the local untracked `bench/format-shape/**/*.vcf.gz`
snapshot present in this worktree at the time of the run.

| Set | n | Faster | Slower | Neutral | Before total | Simplified total | Delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| All inputs | 27 | 18 | 3 | 6 | 48.5999 | 47.0864 | -3.11% |
| Inputs with mean real >= 0.1s | 19 | 17 | 1 | 1 | 48.5033 | 46.9831 | -3.13% |
| FORMAT-hit inputs with mean real >= 0.1s | 15 | 14 | 0 | 1 | 41.5234 | 40.0798 | -3.48% |

Representative FORMAT-heavy deltas:

| Input | Before mean real | Simplified mean real | Delta |
| --- | ---: | ---: | ---: |
| `1000g_chr22_full_genotypes` | 10.9167 | 10.5067 | -3.76% |
| `HG002_CHM13v2.0_v5.0q_smvar` | 1.5633 | 1.3733 | -12.15% |
| `HG002_GRCh38_1_22_v4.2.1_benchmark` | 3.7733 | 3.5767 | -5.21% |
| `HG002_GRCh38_v5.0q_stvar` | 1.8800 | 1.7533 | -6.74% |
| `ccdg_chr22_10k` | 2.5467 | 2.4400 | -4.19% |
| `large_ccdg_likelihood_2048s` | 4.1833 | 4.1000 | -1.99% |
| `large_multiallelic_likelihood_2048s` | 3.0433 | 3.0133 | -0.99% |
| `large_reordered_likelihood_2048s` | 2.7767 | 2.6733 | -3.72% |

Interpretation: the removed width-specific integer vector specializations did
not show a performance benefit on this corpus.  The simpler executor was
slightly faster overall and easier to justify for review.

## Rejected or Deferred Approaches

Width-specific integer parsers were removed because their code-size and review
cost were not supported by the benchmark evidence.

Runtime JIT/code generation remains deferred.  It would avoid some source-level
specialization, but it would add platform, compiler, security, packaging, and
debugging complexity that is hard to justify before the non-JIT planner is
accepted and before HTSlib has a more controlled internal API boundary for VCF
record representation.
