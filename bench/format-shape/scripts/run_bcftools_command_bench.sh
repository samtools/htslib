#!/bin/sh
set -eu

# Broader production-style command benchmark for the VCF FORMAT planner.
#
# The conversion benchmark in run_bcftools_bench.sh measures one important
# path: VCF text -> BCF output via `bcftools view`.  This script intentionally
# exercises a wider set of bcftools command shapes so we can see which workflows
# actually expose FORMAT parse cost:
#
#   view_bcf      full VCF -> BCF conversion
#   view_sites    VCF -> BCF after dropping genotypes with -G
#   query_sites   fixed-column/INFO-oriented query
#   query_format  FORMAT accessor query for a small sample subset
#   stats         bcftools stats
#   filter_gt     FORMAT expression filtering for a small sample subset
#   merge_self    bcftools merge of the input with itself using --force-samples
#
# Each command is run twice, once with HTS_VCF_FORMAT_PLAN=0 and once with
# HTS_VCF_FORMAT_PLAN=1.  Outputs are compared with cmp whenever the command is
# applicable.  FORMAT commands are skipped for sites-only inputs.
#
# Keep the default THREADS_LIST narrow here.  This harness multiplies inputs by
# commands by planner modes, so exhaustive thread scaling belongs in the
# dedicated threaded runner unless a specific command needs investigation.

bcftools=${BCFTOOLS:-bcftools}
inputs=${1:-bench/format-shape/large/bcftools-command-inputs.tsv}
outdir=${OUTDIR:-bench/format-shape/large/results-bcftools-commands}
keep_outputs=${KEEP_OUTPUTS:-1}
threads_list=${THREADS_LIST:-0}
commands=${COMMANDS:-view_bcf view_sites query_sites query_format stats filter_gt}
query_sample_count=${QUERY_SAMPLE_COUNT:-2}
mkdir -p "$outdir"

timings="$outdir/timings.tsv"
checks="$outdir/checks.tsv"
cmds_out="$outdir/commands.tsv"

printf 'name\tcommand\tthreads\tmode\treal\tuser\tsys\n' > "$timings"
printf 'name\tcommand\tthreads\tcomparison\tstatus\n' > "$checks"
printf 'command\tdescription\n' > "$cmds_out"
printf 'view_bcf\tbcftools view --no-version -Ob -l 0\n' >> "$cmds_out"
printf 'view_sites\tbcftools view --no-version -G -Ob -l 0\n' >> "$cmds_out"
printf 'query_sites\tbcftools query fixed site fields\n' >> "$cmds_out"
printf 'query_format\tbcftools query GT for first QUERY_SAMPLE_COUNT samples\n' >> "$cmds_out"
printf 'stats\tbcftools stats\n' >> "$cmds_out"
printf 'filter_gt\tbcftools view -i GT="alt" for first QUERY_SAMPLE_COUNT samples\n' >> "$cmds_out"
printf 'merge_self\tbcftools merge --no-index --force-samples of the input with itself\n' >> "$cmds_out"

run_one()
{
    mode=$1
    command=$2
    threads=$3
    path=$4
    out=$5
    err=$6
    sample_args=$7
    thread_args=
    plan=0

    if [ "$mode" = plan ]; then
        plan=1
    fi
    if [ "$threads" != 0 ]; then
        thread_args="--threads $threads"
    fi

    case "$command" in
        view_bcf)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -Ob -l 0 $thread_args \
                -o "$out" "$path" 2> "$err"
            ;;
        view_sites)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -G -Ob -l 0 $thread_args \
                -o "$out" "$path" 2> "$err"
            ;;
        query_sites)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
                "$path" > "$out" 2> "$err"
            ;;
        query_format)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" query $sample_args -f '%CHROM\t%POS[\t%GT]\n' \
                "$path" > "$out" 2> "$err"
            ;;
        stats)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" stats "$path" > "$out" 2> "$err"
            ;;
        filter_gt)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -Ob -l 0 $thread_args $sample_args \
                -i 'GT="alt"' -o "$out" "$path" 2> "$err"
            ;;
        merge_self)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" merge --no-index --force-samples --no-version -Ob \
                $thread_args -o "$out" "$path" "$path" 2> "$err"
            ;;
        *)
            printf 'unknown command: %s\n' "$command" >&2
            return 1
            ;;
    esac
}

tail -n +2 "$inputs" | while IFS='	' read -r name path source
do
    samples=$("$bcftools" query -l "$path" | awk -v n="$query_sample_count" '
        NR <= n { if (s) s = s "," $0; else s = $0 }
        END { print s }
    ')
    sample_args=
    if [ -n "$samples" ]; then
        sample_args="-s $samples"
    fi

    for command in $commands
    do
        case "$command" in
            query_format|filter_gt|merge_self)
                if [ -z "$sample_args" ]; then
                    for threads in $threads_list
                    do
                        printf '%s\t%s\t%s\tbaseline_vs_plan\tskipped_no_samples\n' \
                               "$name" "$command" "$threads" >> "$checks"
                    done
                    continue
                fi
                ;;
        esac

        for threads in $threads_list
        do
            base_out="$outdir/$name.$command.t$threads.baseline.out"
            plan_out="$outdir/$name.$command.t$threads.plan.out"

            for mode in baseline plan
            do
                err="$outdir/$name.$command.t$threads.$mode.stderr"
                out="$outdir/$name.$command.t$threads.$mode.out"
                run_one "$mode" "$command" "$threads" "$path" "$out" "$err" "$sample_args"

                awk -v name="$name" -v command="$command" \
                    -v threads="$threads" -v mode="$mode" '
                    /^real / { real=$2 }
                    /^user / { user=$2 }
                    /^sys / { sys=$2 }
                    END {
                        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                               name, command, threads, mode,
                               real+0, user+0, sys+0
                    }
                ' "$err" >> "$timings"
            done

            if cmp "$base_out" "$plan_out" >/dev/null 2>&1; then
                printf '%s\t%s\t%s\tbaseline_vs_plan\tok\n' \
                       "$name" "$command" "$threads" >> "$checks"
            else
                printf '%s\t%s\t%s\tbaseline_vs_plan\tDIFF\n' \
                       "$name" "$command" "$threads" >> "$checks"
            fi
            if [ "$keep_outputs" = 0 ]; then
                rm -f "$base_out" "$plan_out"
            fi
        done
    done
done

printf 'wrote %s, %s, and %s\n' "$timings" "$checks" "$cmds_out"
