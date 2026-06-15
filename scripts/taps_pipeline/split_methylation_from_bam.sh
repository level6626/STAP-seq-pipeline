#!/usr/bin/env bash
set -Eeuo pipefail

# Split an existing TAPS aligned BAM by the METH_CODE value embedded in the read
# name, then make per-code alignment-coverage bedGraph/BigWig files.

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-TAPS_27ac_rep1_S7}"
BAM="${BAM:-results/taps_pipeline/${SAMPLE}/${SAMPLE}.bowtie2.paired.sorted.bam}"
OUTDIR="${OUTDIR:-results/taps_pipeline/${SAMPLE}/methylation_split}"
CHROM_SIZES="${CHROM_SIZES:-../data/hg38/hg38.chrom.sizes}"
MIN_MAPQ="${MIN_MAPQ:-0}"
COVERAGE_EXCLUDE_FLAGS="${COVERAGE_EXCLUDE_FLAGS:-2820}"
INCLUDE_UNKNOWN="${INCLUDE_UNKNOWN:-1}"
REQUIRE_BIGWIG="${REQUIRE_BIGWIG:-1}"

LOGDIR="${OUTDIR}/logs"
mkdir -p "${OUTDIR}" "${LOGDIR}"

trap 'echo "ERROR at line ${LINENO}. See ${LOGDIR} for logs." >&2' ERR

log_msg() {
    printf '[%(%F %T)T] %s\n' -1 "$*" | tee -a "${LOGDIR}/taps_methylation_split.progress.log" >&2
}

run_logged() {
    local name="$1"
    shift
    log_msg "START ${name}"
    {
        printf '[%(%F %T)T] command:' -1
        printf ' %q' "$@"
        printf '\n'
        "$@"
        printf '[%(%F %T)T] done\n' -1
    } >"${LOGDIR}/${name}.log" 2>&1
    log_msg "DONE  ${name}"
}

run_bash_logged() {
    local name="$1"
    local command="$2"
    log_msg "START ${name}"
    {
        printf '[%(%F %T)T] command: %s\n' -1 "${command}"
        bash -o pipefail -c "${command}"
        printf '[%(%F %T)T] done\n' -1
    } >"${LOGDIR}/${name}.log" 2>&1
    log_msg "DONE  ${name}"
}

require_file() {
    local path="$1"
    if [[ ! -s "${path}" ]]; then
        echo "Required file is missing or empty: ${path}" >&2
        exit 1
    fi
}

require_command() {
    local name="$1"
    if ! command -v "${name}" >/dev/null 2>&1; then
        echo "Required command is not available in PATH: ${name}" >&2
        exit 1
    fi
}

require_file "${BAM}"
require_file "${CHROM_SIZES}"
require_command samtools
require_command gawk
require_command sort
require_command bedtools
if [[ "${REQUIRE_BIGWIG}" == "1" || "${REQUIRE_BIGWIG}" == "true" ]]; then
    require_command bedGraphToBigWig
fi

MAP_TSV="${OUTDIR}/methylation_code_outputs.tsv"
cat >"${MAP_TSV}" <<EOF
code	label	bam
TTT	100pct	${OUTDIR}/${SAMPLE}.TTT_100pct.bam
AAA	0pct	${OUTDIR}/${SAMPLE}.AAA_0pct.bam
CAT	60pct	${OUTDIR}/${SAMPLE}.CAT_60pct.bam
AGT	40pct	${OUTDIR}/${SAMPLE}.AGT_40pct.bam
TGA	20pct	${OUTDIR}/${SAMPLE}.TGA_20pct.bam
TAG	10pct	${OUTDIR}/${SAMPLE}.TAG_10pct.bam
CTA	1pct	${OUTDIR}/${SAMPLE}.CTA_1pct.bam
ATG	0p1pct	${OUTDIR}/${SAMPLE}.ATG_0p1pct.bam
EOF

if [[ "${INCLUDE_UNKNOWN}" == "1" || "${INCLUDE_UNKNOWN}" == "true" ]]; then
    printf 'unknown\tunknown\t%s/%s.unknown.bam\n' "${OUTDIR}" "${SAMPLE}" >>"${MAP_TSV}"
fi

SUMMARY_TSV="${OUTDIR}/${SAMPLE}.taps_methylation_split.summary.tsv"

run_bash_logged split_bam_by_methylation_code \
    "samtools view -h '${BAM}' | gawk -v map_tsv='${MAP_TSV}' -v summary_tsv='${SUMMARY_TSV}' -v threads='${THREADS}' '
        BEGIN {
            FS = OFS = \"\t\"
            while ((getline line < map_tsv) > 0) {
                if (line !~ /^code\t/) {
                    split(line, f, \"\t\")
                    code = f[1]
                    label[code] = f[2]
                    bam[code] = f[3]
                    cmd[code] = \"samtools view -@ \" threads \" -b -o \\047\" f[3] \"\\047 -\"
                    order[++n] = code
                }
            }
            close(map_tsv)
        }
        /^@/ {
            for (i = 1; i <= n; i++) {
                print \$0 | cmd[order[i]]
            }
            next
        }
        {
            code = \"unknown\"
            if (match(\$1, /(^|\\|)METH_CODE=([^|]+)/, m)) {
                code = m[2]
            }
            if (!(code in cmd)) {
                code = \"unknown\"
            }
            if (code in cmd) {
                print \$0 | cmd[code]
                records[code]++
                flag = \$2 + 0
                if (!and(flag, 4) && !and(flag, 256) && !and(flag, 2048)) {
                    primary_mapped[code]++
                }
            } else {
                skipped++
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                close(cmd[order[i]])
            }
            print \"code\", \"label\", \"bam\", \"bam_records\", \"primary_mapped_records\" > summary_tsv
            for (i = 1; i <= n; i++) {
                code = order[i]
                print code, label[code], bam[code], records[code] + 0, primary_mapped[code] + 0 >> summary_tsv
            }
            if (skipped > 0) {
                print \"skipped\", \"skipped\", \"NA\", skipped, \"NA\" >> summary_tsv
            }
        }
    '"

tail -n +2 "${MAP_TSV}" | while IFS=$'\t' read -r code label split_bam; do
    suffix="${code}_${label}"
    if [[ "${code}" == "${label}" ]]; then
        suffix="${code}"
    fi
    prefix="${OUTDIR}/${SAMPLE}.${suffix}"
    bedgraph="${prefix}.coverage.bedGraph"
    bw="${prefix}.coverage.bw"

    run_logged "index_${suffix}" samtools index "${split_bam}"
    run_logged "flagstat_${suffix}" samtools flagstat "${split_bam}"

    run_bash_logged "coverage_bedgraph_${suffix}" \
        "samtools view -@ ${THREADS} -b -F ${COVERAGE_EXCLUDE_FLAGS} -q ${MIN_MAPQ} '${split_bam}' | bedtools genomecov -bg -ibam stdin | sort -k1,1 -k2,2n > '${bedgraph}'"

    if [[ -s "${bedgraph}" ]] && command -v bedGraphToBigWig >/dev/null 2>&1; then
        run_logged "coverage_bigwig_${suffix}" bedGraphToBigWig "${bedgraph}" "${CHROM_SIZES}" "${bw}"
    elif [[ "${REQUIRE_BIGWIG}" == "1" || "${REQUIRE_BIGWIG}" == "true" ]]; then
        log_msg "SKIP  coverage_bigwig_${suffix}: ${bedgraph} is empty, so no BigWig was created"
    fi
done

log_msg "TAPS methylation split complete for ${SAMPLE}"
column -t -s $'\t' "${SUMMARY_TSV}" | tee -a "${LOGDIR}/taps_methylation_split.progress.log" >&2
