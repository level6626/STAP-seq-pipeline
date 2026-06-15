# STAP-seq Pipeline

This repo contains the working notes and lightweight processing scripts for
Run202 STAP-seq and TAPS-seq data. The STAP discovery scripts still only
discover files whose basename starts with `STAP_`; TAPS processing lives under
`scripts/taps_pipeline`.

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

The STAR 2.7.11b-compatible hg38 index exists at:

```bash
/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24
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
STAR_INDEX_DIR=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24 \
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

## TAPS Methylation Pipeline

TAPS-seq is used here to estimate the real methylation level of each plasmid
barcode. TAPS converts modified cytosines, 5mC/5hmC, to bases read as `T`,
while unmodified `C` remains `C`. For a CpG, this means methylation is observed
as `C>T` on the plus-strand CpG cytosine, or as `G>A` when the read covers the
opposite strand.

The pipeline carries the 17-bp R2 plasmid index into the alignment read name,
aligns the assayed TAPS read, then counts CpG-level conversion grouped by the
first 3 bp methylation code.

Current interpretation:

- `R2[:3]` after orientation is the expected methylation code:
  `TTT` 100%, `AAA` 0%, `CAT` 60%, `AGT` 40%, `TGA` 20%, `TAG` 10%,
  `CTA` 1%, and `ATG` 0.1%.
- `R2[3:17]` is the plasmid/random barcode.
- The emitted read name includes `METH_CODE=`, `METH_EXPECTED=`, `R2=`, and a
  molecule `UMI=` made from the R1 UMI plus the oriented R2 barcode.
- The CpG counter reports both plus-strand CpG `C>T` evidence and opposite
  strand CpG `G>A` evidence.

### TAPS Processing Steps

The wrapper script is `scripts/taps_pipeline/run_taps_pipeline.sh`.

1. Validate inputs:
   - `R1`, `R2`, `R3` are discovered from `RUN_DIR` and `SAMPLE`.
   - `REFERENCE_FASTA` must exist and is indexed with `samtools faidx` if needed.
   - The default STAR index is the STAR `2.7.11b`-compatible hg38 index:
     `/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24`.
2. Parse and tag R2:
   - `scripts/taps_pipeline/prepare_taps_fastqs.py` reads R1/R2/R3 triplets.
   - It orients R2 using `R2_ORIENTATION`, default `forward`.
   - It parses `R2[:3]` as `METH_CODE` and keeps the full `R2[:17]` barcode.
   - It writes tagged R1/R3 FASTQs where the read name contains:
     `METH_CODE=`, `METH_LABEL=`, `METH_EXPECTED=`, `R2=`, and `UMI=`.
   - It writes `<sample>.prepare_taps.stats.tsv`; use this to verify that most
     reads have valid methylation codes and that the chosen R2 orientation is
     correct.
3. Trim reads:
   - `fastp` trims adapters/low-quality sequence on the tagged R1/R3 pair.
   - The default `FASTP_EXTRA=--dont_eval_duplication` avoids high memory use
     from optional fastp duplication-rate estimation.
4. Align:
   - Default genome-wide mode is `ALIGNER=star` and `ALIGN_READS=r3`.
   - `R3` mode matches the older lab mapping branch and worked well for Run202
     TAPS. The old STAR logs mapped `TAPS_27ac_rep1_S7` R3 at roughly 89%
     unique mapping.
   - `ALIGNER=bowtie2` is available for compact reporter/plasmid references.
     Use this only when the reference FASTA really contains the assayed insert
     sequence; a generic oligo workbook FASTA is useful for smoke testing but
     is not expected to capture all genome-wide TAPS reads.
   - `ALIGNER=bismark` is available for bisulfite-style alignment against a
     Bismark genome folder. It supports `ALIGN_READS=r3` and
     `ALIGN_READS=paired`, then emits the same coordinate-sorted BAM interface
     used by the CpG counter.
5. Count CpG conversion:
   - `scripts/taps_pipeline/count_taps_cpg_conversion.py` reads the tagged BAM
     with `pysam`.
   - At every aligned CpG position in `REFERENCE_FASTA`, it counts:
     `converted`, `unconverted`, and `other`.
   - Plus-strand CpG cytosine evidence is `C` unconverted and `T` converted.
   - Opposite-strand CpG evidence is `G` unconverted and `A` converted.
   - Counts are grouped by `METH_CODE`, genomic CpG coordinate, and observed
     strand. The summary collapses over CpGs within each methylation code.

### Run TAPS

Build an optional compact reporter FASTA from the oligo workbook:

```bash
cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin/python \
  scripts/taps_pipeline/write_reporter_fasta_from_oligos.py \
  --oligo-xlsx ../data/meta/STAP_Seq_oligos.xlsx \
  --out results/taps_pipeline/reference/STAP_Seq_oligos.reporter.fa
```

Smoke test the example sample on 2,000 reads:

```bash
THREADS=4 \
MAX_READS=2000 \
SAMPLE=TAPS_27ac_rep1_S7 \
OUTDIR=results/taps_pipeline_smoke/TAPS_27ac_rep1_S7 \
ALIGNER=bowtie2 \
ALIGN_READS=r3 \
REFERENCE_FASTA=results/taps_pipeline/reference/STAP_Seq_oligos.reporter.fa \
R2_ORIENTATION=forward \
scripts/taps_pipeline/run_taps_pipeline.sh
```

For the genome-wide 27ac TAPS run, use STAR on R3 to match the older lab
mapping branch:

```bash
THREADS=12 \
MAX_READS=0 \
SAMPLE=TAPS_27ac_rep1_S7 \
OUTDIR=results/taps_pipeline/TAPS_27ac_rep1_S7 \
ALIGNER=star \
ALIGN_READS=r3 \
REFERENCE_FASTA=../data/hg38/hg38.fa \
STAR_INDEX_DIR=/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24 \
R2_ORIENTATION=forward \
scripts/taps_pipeline/run_taps_pipeline.sh
```

To align the same TAPS reads with Bismark instead, point `REFERENCE_FASTA` at
the FASTA used to create the Bismark genome folder. If
`${BISMARK_GENOME_DIR}/Bisulfite_Genome` is missing and
`BUILD_BISMARK_INDEX=auto` or `1`, the wrapper runs
`bismark_genome_preparation` first.

Single-end R3 Bismark template:

```bash
THREADS=12 \
MAX_READS=0 \
SAMPLE=TAPS_27ac_rep1_S7 \
OUTDIR=results/taps_pipeline/TAPS_27ac_rep1_S7_bismark_r3 \
ALIGNER=bismark \
ALIGN_READS=r3 \
REFERENCE_FASTA=../data/hg38/hg38.fa \
BISMARK_GENOME_DIR=../data/hg38 \
BUILD_BISMARK_INDEX=auto \
BISMARK_PARALLEL=2 \
R2_ORIENTATION=forward \
scripts/taps_pipeline/run_taps_pipeline.sh
```

Paired-end Bismark template:

```bash
THREADS=12 \
MAX_READS=0 \
SAMPLE=TAPS_27ac_rep1_S7 \
OUTDIR=results/taps_pipeline/TAPS_27ac_rep1_S7_bismark_paired \
ALIGNER=bismark \
ALIGN_READS=paired \
REFERENCE_FASTA=../data/hg38/hg38.fa \
BISMARK_GENOME_DIR=../data/hg38 \
BUILD_BISMARK_INDEX=0 \
BISMARK_PARALLEL=2 \
TRIM_R1_UMI=1 \
R2_ORIENTATION=forward \
scripts/taps_pipeline/run_taps_pipeline.sh
```

The wrapper expects `bismark` and `bismark_genome_preparation` in
`${CONDA_ENV_BIN}` or in paths supplied with `BISMARK=` and
`BISMARK_GENOME_PREPARATION=`. The ordinary Bowtie2 index prefix is not enough
for Bismark; the genome folder must contain Bismark's `Bisulfite_Genome`
directory. Use `BISMARK_EXTRA` for chemistry/library-specific flags such as
`--non_directional`.

For libraries where the old notebook removed R2 reads containing the plasmid
motif `TCGGCCTATCATCTGGG`, run with:

```bash
FILTER_R2_MOTIF=TCGGCCTATCATCTGGG scripts/taps_pipeline/run_taps_pipeline.sh
```

### TAPS Output Files

Main outputs:

- `<sample>.prepare_taps.stats.tsv`: raw R2 parsing, methylation-code counts,
  chosen R2 orientation counts, and read discard counts.
- `<sample>.<aligner>.<mode>.sorted.bam`: coordinate-sorted tagged alignment.
- `<sample>.taps_cpg_sites.tsv`: CpG-site conversion table.
- `<sample>.taps_meth_code_summary.tsv`: expected vs observed conversion
  summary for methylation-treatment QC.

Important log files:

- `logs/star_align.log`: exact STAR command and fatal alignment/index errors.
- `logs/<sample>.STAR.Log.final.out`: STAR mapping summary.
- `logs/bismark_align.log`: exact Bismark command and alignment/index errors.
- `logs/fastp_trim.log`: fastp command and trimming status.
- `logs/count_taps_cpg_conversion.log`: CpG-count command and completion status.

### TAPS QC

First check R2 parsing. The expected result is that most reads have one of the
eight designed methylation codes and `r2_orientation_forward` dominates for
Run202:

```bash
cat results/taps_pipeline/TAPS_27ac_rep1_S7/TAPS_27ac_rep1_S7.prepare_taps.stats.tsv
```

If many reads are under `discard_r2_unknown_meth_code`, rerun a small test with:

```bash
R2_ORIENTATION=reverse-complement MAX_READS=200000 scripts/taps_pipeline/run_taps_pipeline.sh
```

or:

```bash
R2_ORIENTATION=both MAX_READS=200000 scripts/taps_pipeline/run_taps_pipeline.sh
```

Then compare the orientation and unknown-code counts.

Next check alignment:

```bash
cat results/taps_pipeline/TAPS_27ac_rep1_S7/logs/TAPS_27ac_rep1_S7.STAR.Log.final.out
cat results/taps_pipeline/TAPS_27ac_rep1_S7/logs/flagstat_sorted.log
```

The default TAPS STAR run keeps unique alignments only
(`--outFilterMultimapNmax 1`), so it cannot be used to estimate how many reads
would have mapped to multiple loci. To audit whether methylation code affects
unique versus multi mapping, rerun STAR while retaining multimappers:

```bash
cd STAP-seq-pipeline
THREADS=16 \
MAX_MULTIMAP=100 \
SOURCE_OUTDIR=results/taps_pipeline/TAPS_27ac_rep1_S7 \
OUTDIR=results/taps_pipeline/TAPS_27ac_rep1_S7_multimap_audit \
scripts/taps_pipeline/run_taps_multimap_audit.sh
```

This reuses the existing tagged and trimmed R3 FASTQ from the full TAPS run,
then writes:

- `TAPS_27ac_rep1_S7.star.r3.multimap100.sorted.bam`: STAR alignment with
  multimappers retained.
- `TAPS_27ac_rep1_S7.mapping_by_meth_code.multimap100.tsv`: per-code mapping
  table. Compare `AAA` and `TTT` using `unique_fraction_vs_input` and
  `multi_fraction_among_bam_mapped`.
- `logs/TAPS_27ac_rep1_S7.STAR.multimap.Log.final.out`: global STAR unique,
  multiple-loci, and too-many-loci rates.

For `TAPS_27ac_rep1_S7`, a very low STAR mapping rate usually means the wrong
STAR index/reference is being used, or an old incompatible STAR index was used.
STAR `2.7.11b` must use:

```bash
/gpfs/data/zhou-lab/yczhang/methylation/data/hg38/STAR_index_2.7.11b_gencode_v24
```

Finally check methylation calibration. `AAA` should be low conversion and `TTT`
should be high conversion. Intermediate codes should generally be ordered by
expected methylation, though local sequence context and coverage can make any
single CpG noisy:

```bash
column -t results/taps_pipeline/TAPS_27ac_rep1_S7/TAPS_27ac_rep1_S7.taps_meth_code_summary.tsv | sed -n '1,12p'
```

Useful summary columns:

- `expected_conversion`: expected in vitro methylation level from the 3-bp code.
- `converted_count`: CpG observations read as TAPS-converted.
- `unconverted_count`: CpG observations read as unconverted.
- `callable_count`: `converted_count + unconverted_count`.
- `conversion_rate`: observed TAPS conversion rate.

For site-level QC, sort CpGs by coverage or inspect a specific methylation code:

```bash
sites=results/taps_pipeline/TAPS_27ac_rep1_S7/TAPS_27ac_rep1_S7.taps_cpg_sites.tsv
head -n 1 "${sites}"
awk 'BEGIN{FS=OFS="\t"} NR>1 && $1=="TTT"' "${sites}" \
  | sort -t $'\t' -k11,11nr \
  | sed -n '1,20p'
```

To check non-CpG background conversion, count TAPS-like conversion at CpH
contexts (`CpA`, `CpC`, and `CpT`). Plus-strand `C>T` and opposite-strand
`G>A` are counted as converted:

```bash
cd STAP-seq-pipeline
scripts/taps_pipeline/run_taps_noncpg_conversion.sh
```

The output is:

```bash
results/taps_pipeline/TAPS_27ac_rep1_S7/TAPS_27ac_rep1_S7.taps_noncpg_context_summary.tsv
```

For the main non-CpG rate by methylation code, inspect the combined `CpH`
rows:

```bash
awk 'BEGIN{FS=OFS="\t"} NR==1 || ($4=="CpH" && $5=="both")' \
  results/taps_pipeline/TAPS_27ac_rep1_S7/TAPS_27ac_rep1_S7.taps_noncpg_context_summary.tsv
```

To also include CpG rows as a direct baseline in the same table:

```bash
INCLUDE_CPG=1 scripts/taps_pipeline/run_taps_noncpg_conversion.sh
```

Interpretation:

- Good R2 QC plus good STAR mapping but poor `TTT`/`AAA` separation suggests a
  methylation/TAPS chemistry issue or a reference/read-orientation issue in the
  CpG counter.
- Poor R2 QC suggests the R2 orientation or read structure is wrong.
- Poor mapping suggests the wrong alignment reference/index, adapter problems,
  or that the assayed TAPS read is not represented by the selected reference.
- Low STAP/TAPS barcode overlap after alignment should be compared with raw R2
  overlap to separate biological/library overlap from alignment losses.

## STAP/TAPS Barcode Overlap

For matched STAP and TAPS libraries, compare overlap using the 17-bp R2 plasmid
barcode, not the full STAP UMI suffix. In the standard STAP BAM, the final
25-bp query-name suffix is `R1_UMI(8) + R2(17)`. In the TAPS BAM, the same
plasmid barcode is stored as `R2=<17bp>`.

Full barcode-overlap run:

```bash
/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin/python \
  scripts/compare_stap_taps_barcodes.py \
  --stap-bam results/standard_tools/STAP_TSS_27ac_rep1_S3/STAP_TSS_27ac_rep1_S3.star.dedup.bam \
  --taps-bam results/taps_pipeline/TAPS_27ac_rep1_S7/TAPS_27ac_rep1_S7.star.r3.sorted.bam \
  --outdir results/barcode_overlap/STAP_TSS_27ac_rep1_S3__TAPS_27ac_rep1_S7 \
  --min-mapq 0
```

Outputs:

- `barcode_counts_stap.tsv`: STAP aligned-record and molecule counts per R2.
- `barcode_counts_taps.tsv`: TAPS aligned-record and molecule counts per R2.
- `barcode_overlap.tsv`: union table with STAP/TAPS presence and counts.
- `barcode_overlap_summary.tsv`: overall and methylation-code-stratified overlap.

To separate barcode-library overlap from alignment effects, compare raw R2 files
directly:

```bash
/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin/python \
  scripts/compare_raw_r2_barcodes.py \
  --stap-r2 ../data/Run202/STAP_TSS_27ac_rep1_S3_R2_001.fastq.gz \
  --taps-r2 ../data/Run202/TAPS_27ac_rep1_S7_R2_001.fastq.gz \
  --outdir results/raw_barcode_overlap/STAP_TSS_27ac_rep1_S3__TAPS_27ac_rep1_S7
```

Raw outputs mirror the alignment overlap outputs:

- `raw_barcode_counts_stap.tsv`
- `raw_barcode_counts_taps.tsv`
- `raw_barcode_overlap.tsv`
- `raw_barcode_overlap_summary.tsv`

Add `--require-known-code` if you want to discard raw R2 reads whose first
3 bp are not one of the eight designed methylation codes before calculating
overlap.

## STAP/TAPS Read Association

After confirming barcode overlap, associate aligned STAP RNA reads with aligned
TAPS DNA reads that share the same 17-bp R2 barcode. This is the same join key
used above: `methylation_code(3 bp) + plasmid_barcode(14 bp)`.

Small capped smoke test:

```bash
/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin/python \
  scripts/associate_stap_taps_reads.py \
  --stap-bam results/standard_tools/STAP_TSS_27ac_rep1_S3/STAP_TSS_27ac_rep1_S3.star.dedup.bam \
  --taps-bam results/taps_pipeline/TAPS_27ac_rep1_S7_paired_bowtie2/TAPS_27ac_rep1_S7.bowtie2.paired.sorted.bam \
  --outdir results/stap_taps_association_smoke/STAP_TSS_27ac_rep1_S3__TAPS_27ac_rep1_S7_paired_bowtie2 \
  --max-stap-records 600000 \
  --max-taps-records 450000 \
  --max-associations 1000
```

Full run template:

```bash
/gpfs/data/zhou-lab/yczhang/miniforge3/envs/stap-standard-tools/bin/python \
  scripts/associate_stap_taps_reads.py \
  --stap-bam results/standard_tools/STAP_TSS_27ac_rep1_S3/STAP_TSS_27ac_rep1_S3.star.dedup.bam \
  --taps-bam results/taps_pipeline/TAPS_27ac_rep1_S7_paired_bowtie2/TAPS_27ac_rep1_S7.bowtie2.paired.sorted.bam \
  --outdir results/stap_taps_association/STAP_TSS_27ac_rep1_S3__TAPS_27ac_rep1_S7_paired_bowtie2
```

Default behavior:

- Uses primary, non-duplicate, mapped alignments with `MAPQ >= 0`.
- Associates `read1` records only, to avoid double-counting paired mates.
- Keeps only the eight known 3-bp methylation codes unless
  `--include-unknown-meth-code` is set.
- Uses a temporary SQLite database in the output directory and removes it when
  complete. Add `--keep-sqlite` to retain it for debugging.

Outputs:

- `stap_taps_read_associations.tsv.gz`: STAP/TAPS read-pair associations.
- `barcode_association_summary.tsv`: per-barcode STAP counts, TAPS counts, and
  possible association-row counts.
- `associate_stap_taps_reads.summary.tsv`: parser/filter metrics and run
  parameters.

## EM-seq Methylation Pipeline

Run188 custom EM-seq triplets are handled by
`scripts/emseq_pipeline/run_emseq_pipeline.sh`. The preprocessing script streams
the synchronized `R1/R2/R3` FASTQs in binary mode, trims the 8-bp R1 UMI from
the emitted R1 by default, and writes paired `R1/R3` FASTQs for Bismark. The R2
read is stored in the read name as metadata:

- `R1[:8]`: UMI.
- `R1[8:]`: genomic/TSS sequence aligned as read 1.
- `R2[:2]`: EM methylation spike-in code:
  `TT=0_pct`, `AA=100_pct`, `GG=40_pct`, `CC=10_pct`, `AT=1_pct`,
  `TA=0.1_pct`; all other 2-bp codes are kept as `Unknown`.
- `R2[2:]`: random plasmid barcode.
- `R3`: reverse genomic mate aligned as read 2.

The emitted read-name format keeps tags parseable and places the molecule UMI
after the final underscore for `umi_tools dedup --extract-umi-method=read_id`:

```text
@Original|METH_CODE=CC|METH_LABEL=10_pct|R2=CCGACTGCTGCGTT|R2_BARCODE=GACTGCTGCGTT|UMI=TNGGTAACGACTGCTGCGTT|UMITOOLS_TNGGTAACGACTGCTGCGTT
```

By default, the molecule UMI is `R1_UMI + R2[2:]`. Override with
`UMI_SOURCE=r1` or `UMI_SOURCE=r1+r2` if the deduplication key should be
narrower or include the methylation code.

Smoke test preprocessing and trimming only, useful before Bismark is installed:

```bash
source /gpfs/data/zhou-lab/yczhang/miniforge3/etc/profile.d/conda.sh
conda activate stap-standard-tools
cd /gpfs/data/zhou-lab/yczhang/methylation/STAP-seq-pipeline

THREADS=2 \
MAX_READS=2000 \
RUN_BISMARK=0 \
SAMPLE=EM_pSTAP_cell_27ac_S14 \
OUTDIR=results/emseq_pipeline_smoke/EM_pSTAP_cell_27ac_S14 \
scripts/emseq_pipeline/run_emseq_pipeline.sh
```

Full hg38 run template:

```bash
THREADS=12 \
BISMARK_PARALLEL=2 \
EXTRACTOR_PARALLEL=4 \
MAX_READS=0 \
SAMPLE=EM_pSTAP_cell_27ac_S14 \
OUTDIR=results/emseq_pipeline/EM_pSTAP_cell_27ac_S14 \
REFERENCE_FASTA=../data/hg38/hg38.fa \
scripts/emseq_pipeline/run_emseq_pipeline.sh
```

The wrapper expects `bismark`, `bismark_genome_preparation`, and
`bismark_methylation_extractor` in `${CONDA_ENV_BIN}` or in paths supplied with
`BISMARK=`, `BISMARK_GENOME_PREPARATION=`, and
`BISMARK_METHYLATION_EXTRACTOR=`. If
`${BISMARK_GENOME_DIR}/Bisulfite_Genome` is missing and
`BUILD_BISMARK_INDEX=auto` or `1`, it runs `bismark_genome_preparation` first.
The existing standard Bowtie2 hg38 index in `../data/hg38/bowtie2_index` is not
the same as a Bismark bisulfite genome index. `THREADS` is used for `fastp` and
`samtools`; set `BISMARK_PARALLEL` and `EXTRACTOR_PARALLEL` separately because
Bismark parallel jobs can consume multiple CPU threads each.

If a run already completed preprocessing and `fastp`, resume from the trimmed
FASTQs with:

```bash
REUSE_TRIMMED_FASTQS=1 \
BUILD_BISMARK_INDEX=0 \
scripts/emseq_pipeline/run_emseq_pipeline.sh
```

For Run188 EM samples whose name contains `500`, the wrapper automatically uses
the compact oligo reference:

```text
/gpfs/data/zhou-lab/dcai/059_DT/1_meCpG_STAP/data_500/data_500.fa
```

For oligo500 runs, the wrapper also creates a right-padded local copy of the
FASTA under `${OUTDIR}/reference/right_padded_oligo` and builds the Bismark
index there. The original 250-bp oligo records end exactly where many R3 reads
align, because the terminal 10 bp contain the oligo barcode. Bismark can report
these as unique Bowtie2 hits and then discard them with
`genomic sequence could not be extracted`; right-padding gives Bismark enough
downstream context while preserving original 1-based oligo coordinates. Disable
this with `PAD_OLIGO_REFERENCE=0` or adjust with `OLIGO_REFERENCE_RIGHT_PAD=`.

EM-seq QC can be run after alignment with:

```bash
SAMPLE=EM_pSTAP_cell_500_S13 \
OUTDIR=results/emseq_pipeline/EM_pSTAP_cell_500_S13_padded \
scripts/emseq_pipeline/run_emseq_qc.sh
```

For oligo500 samples, QC also checks whether the mapped oligo's terminal barcode
matches the R3/read2 barcode. The default barcode length is 10 bp and can be
changed with `OLIGO_BARCODE_LENGTH=`. The non-CpG conversion table includes CpG
rows as a baseline by default.
