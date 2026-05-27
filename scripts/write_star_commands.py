#!/usr/bin/env python3
"""Generate STAR commands for genomic STAP TSS R3 reads."""

import argparse
import csv
from pathlib import Path


DEFAULT_STAR_INDEX = "/gpfs/data/zhou-lab/dcai/data/hg38/STAR_index/STAR"
DEFAULT_GTF = "/gpfs/data/zhou-lab/dcai/data/hg38/gencode.v24.annotation.gtf"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--output-root", default="results/genomic_star")
    parser.add_argument("--star-index", default=DEFAULT_STAR_INDEX)
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--threads", type=int, default=6)
    args = parser.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    output_root = Path(args.output_root)

    rows = []
    with open(args.manifest) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["assay_branch"] == "tss_genomic":
                rows.append(row)

    with out.open("w") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n\n")
        handle.write('command -v STAR >/dev/null || { echo "STAR is not on PATH" >&2; exit 1; }\n')
        handle.write(f'mkdir -p "{output_root}"\n\n')
        for row in rows:
            sample = row["sample"]
            r3 = row["R3"]
            sample_dir = output_root / sample
            prefix = sample_dir / sample
            handle.write(f'mkdir -p "{sample_dir}"\n')
            handle.write(
                "STAR "
                "--runMode alignReads "
                "--genomeLoad NoSharedMemory "
                "--readFilesCommand zcat "
                "--outSAMstrandField intronMotif "
                f"--runThreadN {args.threads} "
                f'--genomeDir "{args.star_index}" '
                f'--readFilesIn "{r3}" '
                "--clip5pNbases 3 "
                f'--sjdbGTFfile "{args.gtf}" '
                f'--outFileNamePrefix "{prefix}" '
                "--outSAMreadID Number "
                "--outSAMtype BAM SortedByCoordinate "
                "--outFilterMultimapNmax 1 "
                "--outFilterMismatchNmax 10 "
                "--outFilterMismatchNoverLmax 0.1 "
                "--outFilterMatchNminOverLread 0.3 "
                "--outFilterScoreMinOverLread 0.3\n\n"
            )


if __name__ == "__main__":
    main()
