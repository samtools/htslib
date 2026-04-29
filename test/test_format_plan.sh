#!/bin/sh
set -eu

test_view=${TEST_VIEW:-./test/test_view}
inputs=${1:-"test/format-plan-edge.vcf test/format-plan-header-mismatch.vcf test/format-plan-composable.vcf test/format-plan-gt-header-shape.vcf test/format-plan-cache.vcf test/format-plan-profitability.vcf"}
tmpdir=${TMPDIR:-/tmp}
base=${tmpdir}/hts-format-plan-base.$$
plan=${tmpdir}/hts-format-plan-plan.$$
interp=${tmpdir}/hts-format-plan-interp.$$
stats=${tmpdir}/hts-format-plan-stats.$$
interp_stats=${tmpdir}/hts-format-plan-interp-stats.$$

trap 'rm -f "$base" "$plan" "$interp" "$stats" "$interp_stats"' EXIT HUP INT TERM

for input in $inputs
do
    env HTS_VCF_FORMAT_PLAN=0 "$test_view" -b -l 0 "$input" > "$base"
    env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 "$test_view" -b -l 0 "$input" > "$plan" 2> "$stats"
    env HTS_VCF_FORMAT_PLAN=interp HTS_VCF_FORMAT_PLAN_STATS=1 "$test_view" -b -l 0 "$input" > "$interp" 2> "$interp_stats"
    cmp "$base" "$plan"
    cmp "$base" "$interp"
    case "$input" in
        *format-plan-cache.vcf)
            grep -q 'attempts=21 hits=21 fallback=0 ' "$stats"
            grep -q 'attempts=21 hits=21 fallback=0 ' "$interp_stats"
            ;;
        *format-plan-profitability.vcf)
            grep -q 'attempts=1 hits=0 fallback=1 ' "$stats"
            grep -q 'attempts=1 hits=0 fallback=1 ' "$interp_stats"
            ;;
    esac
    cat "$stats"
    cat "$interp_stats"
done

for samples in S1,S3 S2 ^S2
do
    for input in test/format-plan-composable.vcf test/format-plan-edge.vcf
    do
        env HTS_VCF_FORMAT_PLAN=0 "$test_view" -b -l 0 -s "$samples" "$input" > "$base"
        env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 "$test_view" -b -l 0 -s "$samples" "$input" > "$plan" 2> "$stats"
        env HTS_VCF_FORMAT_PLAN=interp HTS_VCF_FORMAT_PLAN_STATS=1 "$test_view" -b -l 0 -s "$samples" "$input" > "$interp" 2> "$interp_stats"
        cmp "$base" "$plan"
        cmp "$base" "$interp"
        cat "$stats"
        cat "$interp_stats"
    done
done
