# STAP-seq Pipeline

This repo contains the working notes and lightweight processing scripts for
Run202 STAP-seq data. Run202 also contains TAPS-seq files; those are deliberately
excluded here. Every script discovers only files whose basename starts with
`STAP_`.

## Run202 Layout

Raw data are available through the symlink:

```bash
../data/Run202 -> /gpfs/data/zhou-lab/dcai/data/Run202
```

There are two STAP library layouts in Run202:

1. `STAP_TSS_*`: three FASTQs per sample, `R1/R2/R3`.
2. `STAP_scTSSV2C_*`: two FASTQs per sample, `R1/R2`.

The current scripts focus on `STAP_TSS_*`. The `STAP_scTSSV2C_*` files are
inventoried but not processed by the TSS commands.

## Processing Decision

The `STAP_TSS_*` triplets are not all the same assay branch:

- `STAP_TSS_500_oligos_*` and `STAP_TSS_6_oligos_*` match compact oligo design
  FASTAs from the older lab analysis:
  - `/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa`
  - `/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_6/data_6.fa`
- `STAP_TSS_27ac_*`, `STAP_TSS_BasPromoter_*`, and `STAP_TSS_SE_*` do not match
  those oligo FASTAs by k-mer fingerprinting. They follow the older genomic
  STAP TSS branch, where `R3` is mapped to hg38 with STAR.

For the `STAP_TSS_*` triplets, the current read interpretation is:

- `R1`: 8-bp RNA UMI followed by the TSS-starting sequence.
- `R2`: 17-bp molecule/plasmid index. The first 3 bp encode methylation level,
  followed by a 14-bp random barcode. The known methylation codes are `TTT`
  100%, `AAA` 0%, `CAT` 60%, `AGT` 40%, `TGA` 20%, `TAG` 10%, `CTA` 1%, and
  `ATG` 0.1%.
- `R3`: paired from the opposite side toward R1. The oligo-control libraries
  and `SE` contain the DNA barcode in this read.

The shared hg38 STAR index exists at:

```bash
/gpfs/data/zhou-lab/dcai/data/hg38/STAR_index/STAR
```

The older V3 script used this annotation:

```bash
/gpfs/data/zhou-lab/dcai/data/hg38/gencode.v24.annotation.gtf
```

## Quick Start

Create a STAP-only manifest:

```bash
python3 scripts/discover_stap_run.py \
  --run-dir ../data/Run202 \
  --out results/run202_stap_manifest.tsv
```

Quantify the oligo-control STAP TSS samples using the known design FASTAs:

```bash
python3 scripts/quantify_oligo_tss.py \
  --sample STAP_TSS_500_oligos_S2 \
  --run-dir ../data/Run202 \
  --reference /gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa \
  --outdir results/oligo_tss \
  --r2-orientation forward

python3 scripts/quantify_oligo_tss.py \
  --sample STAP_TSS_6_oligos_S1 \
  --run-dir ../data/Run202 \
  --reference /gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_6/data_6.fa \
  --outdir results/oligo_tss \
  --r2-orientation forward
```

If the methylation code appears reverse-complemented, rerun with
`--r2-orientation reverse-complement` and compare `unknown_methylation_code` in
the per-sample `summary.tsv`.

## Oligo Alignment And Counts

`scripts/quantify_oligo_tss.py` does not use a genome aligner. The oligo
references are short enough that the script uses exact k-mer matching:

1. Load the design FASTA and build a non-ambiguous k-mer index from both strands.
   The default k-mer length is 25 bp.
2. For each read triplet, scan `R1` first, then `R3`, for the first exact k-mer
   present in the design index.
3. Convert the matching read offset and reference k-mer offset into a TSS
   position on the oligo.
4. Parse `R1[:8]` as the RNA UMI.
5. Parse oriented `R2` as:
   - `R2[:3]`: methylation code
   - `R2[3:]`: 14-bp random plasmid/molecule barcode
6. Count reads by `sample`, `oligo_id`, `position`, `strand`, and
   `methylation_code`.

Output columns in `tss_counts.tsv`:

- `raw_count`: number of assigned read triplets in that bin.
- `dedup_count`: number of unique molecule keys in that bin, where the molecule
  key is `R1_UMI:R2_oriented`. This collapses reads with the same RNA UMI and
  same 17-bp R2 molecule index.

The script writes `summary.tsv` and `progress.tsv`. Per-read assignments are not
written by default because full runs are large; add `--write-assignments` only
when debugging individual reads.

Generate, but do not run, STAR commands for genomic STAP TSS samples:

```bash
python3 scripts/write_star_commands.py \
  --manifest results/run202_stap_manifest.tsv \
  --out scripts/run202_stap_genomic_star.sh
```

Then run the generated shell script on a node where `STAR` is available.

## Standard-Tools Comparator Pipeline

`scripts/standard_tools/run_stap_standard_tools.sh` is an independent
comparator workflow built around common bioinformatics tools. It is intended for
STAP TSS triplets only, not `TAPS_*` files.

Workflow:

1. Subset reads when `MAX_READS` is set, for smoke testing.
2. Run `umi_tools extract` on `R1/R2`.
   - Removes the 8-bp RNA UMI from the 5' end of `R1`.
   - Appends the combined `R1_UMI + R2_17bp_barcode` to the read name.
3. Synchronize `R3` headers to the UMI-tagged `R1` headers.
   `umi_tools extract` is paired-end rather than triplet-aware, so this helper
   keeps `R1/R3` paired-read names identical for Bowtie2/STAR and
   `umi_tools dedup`.
4. Trim adapters and low-quality bases with `fastp`.
   The script does not enable 5' quality trimming, so the post-UMI first base of
   `R1` remains the assayed TSS base.
5. Align `R1/R3` with either:
   - `ALIGNER=bowtie2` for reporter/plasmid/oligo references without splicing.
   - `ALIGNER=star` for genome references where splicing is expected.
6. Coordinate-sort and index the BAM with `samtools`.
7. Deduplicate with `umi_tools dedup --paired`.
8. Extract the 5'-most mapped position of read 1 into BED, then create
   bedGraph and BigWig TSS coverage.
9. Optionally count read1 TSS events over candidate windows with
   `bedtools coverage -counts`.

Small smoke test, using only the first 2,000 `STAP_TSS_500_oligos_S2` triplets:

```bash
srun --jobid=11320722 --pty bash
source /gpfs/data/zhou-lab/yczhang/miniforge3/etc/profile.d/conda.sh
conda activate stap-standard-tools
cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

THREADS=4 \
MAX_READS=2000 \
SAMPLE=STAP_TSS_500_oligos_S2 \
OUTDIR=results/standard_tools_smoke/STAP_TSS_500_oligos_S2 \
ALIGNER=bowtie2 \
REFERENCE_FASTA=/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa \
scripts/standard_tools/run_stap_standard_tools.sh
```

Full Bowtie2 run template:

```bash
THREADS=16 \
MAX_READS=0 \
SAMPLE=STAP_TSS_500_oligos_S2 \
OUTDIR=results/standard_tools/STAP_TSS_500_oligos_S2 \
ALIGNER=bowtie2 \
REFERENCE_FASTA=/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa \
scripts/standard_tools/run_stap_standard_tools.sh
```

STAR/genome run template:

```bash
THREADS=16 \
MAX_READS=0 \
SAMPLE=STAP_TSS_27ac_rep1_S3 \
OUTDIR=results/standard_tools/STAP_TSS_27ac_rep1_S3 \
ALIGNER=star \
STAR_INDEX_DIR=/gpfs/data/zhou-lab/dcai/data/hg38/STAR_index/STAR \
CHROM_SIZES=/path/to/hg38.chrom.sizes \
scripts/standard_tools/run_stap_standard_tools.sh
```

For large genome-wide libraries, `fastp` duplication-rate evaluation can use
more memory than a small interactive allocation provides. The standard-tools
script disables that optional calculation by default with
`FASTP_EXTRA=--dont_eval_duplication`. This does not remove reads; PCR
deduplication still happens later with `umi_tools dedup`.

If a genome run was interrupted after UMI extraction, reuse the completed
intermediate FASTQs instead of rebuilding them:

```bash
REUSE_UMI_FASTQS=1 \
FASTP_EXTRA=--dont_eval_duplication \
scripts/standard_tools/run_stap_standard_tools.sh
```

`run_star.sh` is configured this way for the existing
`STAP_TSS_27ac_rep1_S3` recovery run. Start it on a compute node with enough
memory for the hg38 STAR index and coordinate sorting, for example:

```bash
salloc --mem=96G --cpus-per-task=12
bash run_star.sh
```

Trimmed FASTQs are now written to temporary filenames and renamed only after
`fastp` succeeds, so interrupted output cannot be mistaken for a completed
trimmed FASTQ.

The older shared hg38 STAR index at
`/gpfs/data/zhou-lab/dcai/data/hg38/STAR_index/STAR` was generated with an older
STAR genome format and cannot be loaded by STAR `2.7.11b`. Build a compatible
index once in the writable project data directory:

```bash
salloc --mem=96G --cpus-per-task=12
source /gpfs/data/zhou-lab/yczhang/miniforge3/etc/profile.d/conda.sh
conda activate stap-standard-tools
cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

THREADS=12 scripts/reference/build_hg38_star_index.sh
```

The build script:

1. Downloads UCSC `hg38.fa.gz`.
2. Validates and decompresses the FASTA.
3. Uses the existing shared `gencode.v24.annotation.gtf`.
4. Generates `data/hg38/STAR_index_2.7.11b_gencode_v24`.

`run_star.sh` points at that new index. After the one-time build completes, run:

```bash
bash run_star.sh
```

Candidate-window counting:

```bash
CANDIDATE_WINDOWS=/path/to/windows.bed \
scripts/standard_tools/run_stap_standard_tools.sh
```

Main outputs:

- `${SAMPLE}.${ALIGNER}.sorted.bam`: sorted, indexed alignment before UMI
  deduplication.
- `${SAMPLE}.${ALIGNER}.dedup.bam`: sorted, indexed alignment after
  `umi_tools dedup`.
- `${SAMPLE}.read1_tss.bed`: one deduplicated read1 TSS event per line.
- `${SAMPLE}.read1_tss.bedGraph`: read1 TSS coverage.
- `${SAMPLE}.read1_tss.bw`: BigWig version of the TSS coverage.
- `logs/`: stdout/stderr for every tool step plus a progress log.

By default, `DEDUP_WRITE_STATS=0` because `umi_tools 1.1.5 --output-stats` can
fail with newer pandas releases. The environment file pins `pandas<2`; if your
existing environment still has pandas 2.x, either keep the default or update it:

```bash
conda install -n stap-standard-tools 'pandas<2'
```

### Split Existing BAM By Methylation Code

The merged standard-tools TSS BED/bedGraph/BigWig files combine all R2
methylation codes. For STAP oligo/reporter data, split the existing aligned BAM
afterward by parsing the first 3 bases of the R2 barcode embedded in the read
name.

For the full oligo500 run:

```bash
source /gpfs/data/zhou-lab/yczhang/miniforge3/etc/profile.d/conda.sh
conda activate stap-standard-tools
cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

THREADS=8 \
SAMPLE=STAP_TSS_500_oligos_S2 \
BAM=results/standard_tools/STAP_TSS_500_oligos_S2/STAP_TSS_500_oligos_S2.bowtie2.dedup.bam \
OUTDIR=results/standard_tools/STAP_TSS_500_oligos_S2/methylation_split \
CHROM_SIZES=results/standard_tools/STAP_TSS_500_oligos_S2/reference/chrom.sizes \
scripts/standard_tools/split_methylation_from_bam.sh
```

This creates per-code BAM, BED, bedGraph, and BigWig files for:

- `TTT_100pct`
- `AAA_0pct`
- `CAT_60pct`
- `AGT_40pct`
- `TGA_20pct`
- `TAG_10pct`
- `CTA_1pct`
- `ATG_0p1pct`
- `unknown`

The parser assumes the read name suffix is `R1_UMI(8bp) + R2(17bp)`, which is
what `run_stap_standard_tools.sh` writes. If you split a BAM generated with only
the 17-bp R2 barcode in the read name, set `APPENDED_R1_UMI_LENGTH=0`.

## Barcode-First Oligo Pipeline

The current experiment-side interpretation for oligo/reporter STAP libraries is:

- `R1`: 8-bp RNA UMI, followed immediately by the TSS-starting sequence.
- `R2`: 3-bp methylation code, then 14-bp secondary random UMI.
- `R3`: 5' DNA oligo barcode identifying the plasmid/oligo variant.
- Oligo barcode dictionary:
  `../data/meta/STAP_Seq_oligos.xlsx`.

Use `scripts/barcode_pipeline/run_stap_barcode_pipeline.sh` for this workflow.
It streams `R1/R2/R3` together, discards reads that do not match the methylation
and oligo dictionaries, and writes tagged paired `R1/R3` FASTQs for alignment.

The output read name format is SAM-safe:

```text
@OriginalReadName|METH=0.1%|OLIGO=NativeTSS:6749:CCATGCACAC_TTGTCAGTTATGTTAGGGGATA
```

The final underscore-delimited suffix is the 22-bp combined UMI:
`R1_UMI(8 bp) + R2_secondary_UMI(14 bp)`. This lets `umi_tools dedup --paired`
use the combined UMI natively. The `METH` and `OLIGO` tags are placed before the
final underscore so Bowtie2 preserves them in the BAM query name.

Small smoke test:

```bash
source /gpfs/data/zhou-lab/yczhang/miniforge3/etc/profile.d/conda.sh
conda activate stap-standard-tools
cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

THREADS=4 \
MAX_READS=2000 \
SAMPLE=STAP_TSS_500_oligos_S2 \
OUTDIR=results/barcode_pipeline_smoke/STAP_TSS_500_oligos_S2 \
REFERENCE_FASTA=/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa \
BARCODE_ORIENTATION=both \
BARCODE_SEARCH_BASES=0 \
scripts/barcode_pipeline/run_stap_barcode_pipeline.sh
```

Full oligo500 template:

```bash
THREADS=16 \
MAX_READS=0 \
SAMPLE=STAP_TSS_500_oligos_S2 \
OUTDIR=results/barcode_pipeline/STAP_TSS_500_oligos_S2 \
REFERENCE_FASTA=/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa \
BARCODE_ORIENTATION=both \
BARCODE_SEARCH_BASES=0 \
scripts/barcode_pipeline/run_stap_barcode_pipeline.sh
```

Full oligo6 template:

```bash
THREADS=8 \
MAX_READS=0 \
SAMPLE=STAP_TSS_6_oligos_S1 \
OUTDIR=results/barcode_pipeline/STAP_TSS_6_oligos_S1 \
REFERENCE_FASTA=/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_6/data_6.fa \
BARCODE_ORIENTATION=both \
BARCODE_SEARCH_BASES=0 \
scripts/barcode_pipeline/run_stap_barcode_pipeline.sh
```

Important options:

- `BARCODE_ORIENTATION=both`: test showed oligo500 R3 barcodes match the
  reverse-complement of the Excel barcode at offset 0.
- `MAX_BARCODE_MISMATCHES=1`: allows one barcode mismatch. Ambiguous one-mismatch
  matches are discarded.
- `BARCODE_SEARCH_BASES=0`: require the barcode at the first base of R3. Increase
  this only if the barcode can be shifted downstream.
- `KEEP_R3_BARCODE=1`: default behavior keeps the matched R3 barcode for
  alignment. This mapped better in a 2,000-read oligo500 smoke test against
  `data_500.fa`, which appears to include barcode sequence. Set to `0` only if
  the Bowtie2 reference lacks the barcode sequence.

Main outputs:

- `${SAMPLE}.demux.stats.tsv`: demultiplexing and discard summary.
- `${SAMPLE}.oligo_metadata.tsv`: loaded oligo metadata and generated `Oligo_ID`
  values.
- `${SAMPLE}.bowtie2.sorted.bam`: sorted alignment before UMI deduplication.
- `${SAMPLE}.bowtie2.dedup.bam`: paired UMI-deduplicated BAM.
- `${SAMPLE}.tss_by_oligo_meth.tsv`: final table with columns
  `Oligo_ID`, `Meth_State`, `Chromosome`, `TSS_Position`, and `Count`.

## Important Separation

Do not process `TAPS_*` files with this repo. TAPS-seq has different chemistry
and should get its own manifest and processing workflow.
