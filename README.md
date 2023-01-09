# QRS - Quick Reads Stats

Quick computation of alignment statistic from BAM files.

## Usage

```bash
Usage:
  CSQ Selector [options] [bam ...]

Arguments:
  [bam ...]        input BAM/CRAM file(s). Glob pattern allowed

Options:
  -h, --help
  -o, --out=OUT              Output file (default: QuickRunStats)
  -s, --subset=SUBSET        Read only N reads from each file. Use -1 for all reads (default: 1000000)
  -q, --minq=MINQ            Min base quality for stats reporting (default: 30)
```

## Description

The tool subsets N random reads from input BAM(s) and report essential statistic about the reads as well as organization of the run (sequencer, run, flow-cell, lane) based on Illumina reads name. The fraction of bases with quality above `--minq` is reported as well.

## Output

Sample names are obtained from the header of BAM file(s) directly based on the @RG fields. The tool then output:

- 1 file with alignment stats `<outprefix>.reads_stat.tsv`
- 2 files for sequencing run, one containing the instrument name(s), run id(s), and flowcell id(s) per sample `<outprefix>.rundetails.tsv`; one containing a summary of how many runs, flow cells and lanes where used for each sample `<outprefix>.per_sample_run.tsv`
- 1 file per sample containing che GC fractions for all reads analyzed. These are saved in a folder named `<outprefix>.GCstats`

### read_stats.tsv

Reports statistics collected across analysed reads.

| Column | Description |
|--------|-------------|
| SAMPLE | Sample id from the BAM file |
| MEAN_READLEN | Mean read length across the sampled reads |
| MIN_READLEN | Minimum read length across the samples reads |
| MEAN_BASEQ | Mean base quality across all bases in the sampled reads |
| MEDIAN_BASEQ | Median base quality across all bases in the sampled reads |
| PERC_BASES_ABOVE_Q{minQ} | Fraction of bases above the value set by `--minQ` |
| MEAN_MAPQ | Mean mapping quality across the sampled reads |
| MEDIAN_MAPQ | Median mapping quality across the sampled reads |
| DUP_RATE | Fraction of duplicated reads across sampled reads (based on read flags) |
| MEAN_READ_GC_PERC | Mean GC perc across sampled reads |

### rundetails.tsv

Reports for each sample the corresponding run ids (instrument, run, flowcell, lane number)

| Column | Description |
|--------|-------------|
| SAMPLE | Sample id from the BAM file |
| INSTRUMENT | Instrument id extracted from read names |
| RUN | Run id extracted from read names |
| FLOWCELL | Flowcell id extracted from read names |
| LANE | Lane number extracted from read names |

### per_sample_run.tsv

Describes how each sample sequencing was distributed across runs, flowcells and lanes.

| Column | Description |
|--------|-------------|
| SAMPLE | Sample id from the BAM file |
| N_RUNS | Number of unique sequencing runs detected for the sample based on read names |
| N_FLOWCELLS | Number of unique flowcells detected for the sample based on read names |
| N_LANES | Number of unique lanes detected for the sample based on read names |
