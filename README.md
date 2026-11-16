# QRS - Quick Reads Stats

Quick computation of read-level and alignment statistics from BAM/CRAM and FASTQ files.

## Usage

```bash
Usage:
  qrs [options]

Options:
  -h, --help
  -b, --bam=BAM              Input BAM/CRAM file(s). Can be repeated or use glob patterns
  -1, --fastq1=FASTQ1        Input R1 (or single-end) FASTQ(.gz) file
  -2, --fastq2=FASTQ2        Optional input R2 FASTQ(.gz) file for paired-end reads
  -n, --sample=SAMPLE        Sample name for FASTQ input (default: derived from --fastq1 filename)
  -o, --out=OUT              Output file prefix (default: QuickRunStats)
  -s, --subset=SUBSET        Read only N random reads per sample. Use -1 for all reads (default: 1000000)
  -q, --minq=MINQ            Min base quality for stats reporting (default: 30)
  --write_read_gc            Write per-read GC statistics (disabled by default)
  --write_run_stats          Write run details and per-sample run summary (disabled by default)
```

## Description

`qrs` subsets `N` random reads from input BAM/CRAM file(s) (`--bam`) or FASTQ/FASTQ.gz file(s) (`--fastq1` and optional `--fastq2`) and reports essential read-level and alignment statistics:

- **BAM/CRAM mode (`--bam`)**: Samples `N` reads across canonical chromosomes (`(chr)?[0-9X]+`) via stratified multi-window index jumps weighted by per-chromosome mapped read counts, computing read length, base quality, mapping quality, duplication rate, proper-pair rate, GC content, and CIGAR-based per-read error rate and indel event count. Sample names are extracted from `@RG` (`SM`) header tags.
- **FASTQ mode (`--fastq1` [`--fastq2`])**: Samples `N` random reads uniformly across uncompressed or gzip/bgzip-compressed FASTQ file(s) via 1-pass reservoir sampling (or streams all reads when `-s -1`), computing read length, base quality, GC content, and — when `--fastq2` is provided — the fraction of reads having a mate pair (`PROPER_PAIR_RATE`). Alignment-only metrics are reported as `NA`.

## Output

The tool outputs:

- `<outprefix>.read_stats.tsv`: summary read-level and alignment statistics per sample.
- With `--write_read_gc`: 1 file per sample containing per-read GC fractions streamed to `<outprefix>.GC_stats/<sample>.GC_stats.tsv` (disabled by default; `MEAN_READ_GC_PERC` is always reported in `read_stats.tsv`).
- With `--write_run_stats`: 2 sequencing-run files extracted from Illumina read names — `<outprefix>.rundetails.tsv` (instrument, run, flowcell, lane per sample) and `<outprefix>.per_sample_run.tsv` (counts of unique runs, flowcells, and lanes per sample). Disabled by default.

### read_stats.tsv

Reports statistics collected across analyzed reads.

| Column                   | Description                                                                                                       |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------- |
| SAMPLE                   | Sample ID (from BAM `@RG` `SM` tag, `--sample`, or derived from `--fastq1` filename)                              |
| MEAN_READLEN             | Mean read length across sampled reads                                                                             |
| MIN_READLEN              | Minimum read length across sampled reads                                                                          |
| MEAN_BASEQ               | Mean base quality across all bases in sampled reads                                                               |
| MEDIAN_BASEQ             | Median base quality across all bases in sampled reads                                                             |
| PERC_BASES_ABOVE_Q{minQ} | Fraction of bases with quality >= `--minq`                                                                        |
| MEAN_MAPQ                | Mean mapping quality across sampled reads (`NA` in FASTQ mode)                                                    |
| MEDIAN_MAPQ              | Median mapping quality across sampled reads (`NA` in FASTQ mode)                                                  |
| DUP_RATE                 | Fraction of duplicate reads based on alignment flags (`NA` in FASTQ mode)                                         |
| PROPER_PAIR_RATE         | Fraction of properly paired reads (from `0x2` flag in BAM mode; fraction of reads with a mate in paired FASTQ; `NA` for single-end FASTQ) |
| MEAN_READ_GC_PERC        | Mean per-read GC fraction across sampled reads                                                                    |
| MEAN_ERROR_RATE          | Mean per-read error rate (`edit_distance / aligned_query_bases` from CIGAR and `NM`; `NA` in FASTQ mode)          |
| MEAN_INDELS              | Mean number of indel events (`I` + `D` CIGAR operations) per aligned read (`NA` in FASTQ mode)                    |

### rundetails.tsv (`--write_run_stats`)

Reports for each sample the corresponding run IDs (instrument, run, flowcell, lane number).

| Column     | Description                             |
| ---------- | --------------------------------------- |
| SAMPLE     | Sample ID                               |
| INSTRUMENT | Instrument ID extracted from read names |
| RUN        | Run ID extracted from read names        |
| FLOWCELL   | Flowcell ID extracted from read names   |
| LANE       | Lane number extracted from read names   |

### per_sample_run.tsv (`--write_run_stats`)

Describes how each sample's sequencing was distributed across runs, flowcells, and lanes.

| Column      | Description                                                                  |
| ----------- | ---------------------------------------------------------------------------- |
| SAMPLE      | Sample ID                                                                    |
| N_RUNS      | Number of unique sequencing runs detected for the sample based on read names |
| N_FLOWCELLS | Number of unique flowcells detected for the sample based on read names       |
| N_LANES     | Number of unique lanes detected for the sample based on read names           |

## Docker

Build the executable inside a multi-stage image:

```bash
docker build -t qrs .
docker run --rm qrs --help
```

You can use the prebuilt Docker image from GitHub Container Registry:

```bash
docker run --rm ghcr.io/htgenomeanalysisunit/quick_reads_stats:v0.2.1 --help
```
