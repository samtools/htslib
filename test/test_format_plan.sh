#!/bin/sh
set -eu

test_view=${TEST_VIEW:-./test/test_view}
input=${1:-test/format-plan-edge.vcf}
tmpdir=${TMPDIR:-/tmp}
base=${tmpdir}/hts-format-plan-base.$$
plan=${tmpdir}/hts-format-plan-plan.$$
stats=${tmpdir}/hts-format-plan-stats.$$

trap 'rm -f "$base" "$plan" "$stats"' EXIT HUP INT TERM

env HTS_VCF_FORMAT_PLAN=0 "$test_view" -b -l 0 "$input" > "$base"
env HTS_VCF_FORMAT_PLAN=1 HTS_VCF_FORMAT_PLAN_STATS=1 "$test_view" -b -l 0 "$input" > "$plan" 2> "$stats"
cmp "$base" "$plan"
cat "$stats"
