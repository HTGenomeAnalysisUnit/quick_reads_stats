import argparse
import strformat
import strutils
from ./utils import log

var p = newParser("qrs"):
    help("Compute read-level and alignment statistics from BAM/CRAM or FASTQ files")
    option("-b", "--bam", multiple=true, help="Input BAM/CRAM file(s). Can be repeated or use glob patterns")
    option("-1", "--fastq1", help="Input R1 (or single-end) FASTQ(.gz) file")
    option("-2", "--fastq2", help="Optional input R2 FASTQ(.gz) file for paired-end reads")
    option("-n", "--sample", help="Sample name for FASTQ input (default: derived from --fastq1 filename)")
    option("-o", "--out", help="Output file prefix", default=some("QuickRunStats"))
    option("-s", "--subset", help="Read only N random reads per sample. Use -1 for all reads", default=some("1000000"))
    option("-q", "--minq", help="Min base quality for stats reporting", default=some("30"))
    flag("--write_read_gc", help="Write per-read GC statistics (disabled by default)")
    flag("--write_run_stats", help="Write run details and per-sample run summary (disabled by default)")

proc parseCmdLine*(): ref =
    try:
        result = p.parse()
        let hasBam = result.bam.len > 0
        let hasFq1 = result.fastq1.len > 0
        let hasFq2 = result.fastq2.len > 0
        if not hasBam and not hasFq1:
            raise UsageError.newException("Either --bam or --fastq1 must be specified")
        if hasBam and (hasFq1 or hasFq2):
            raise UsageError.newException("Cannot specify both --bam and --fastq1/--fastq2 at the same time")
        if hasFq2 and not hasFq1:
            raise UsageError.newException("--fastq2 requires --fastq1 to also be specified")

        var subsetVal: int
        try:
            subsetVal = parseInt(result.subset)
        except ValueError:
            raise UsageError.newException(fmt"Invalid integer for --subset: {result.subset}")
        if subsetVal != -1 and subsetVal <= 0:
            raise UsageError.newException("--subset must be -1 or a positive integer")

        var minqVal: int
        try:
            minqVal = parseInt(result.minq)
        except ValueError:
            raise UsageError.newException(fmt"Invalid integer for --minq: {result.minq}")
        if minqVal < 0:
            raise UsageError.newException("--minq must be >= 0")
    except ShortCircuit as e:
        if e.flag == "argparse_help":
            echo p.help
            quit QuitSuccess
    except UsageError:
        stderr.writeLine getCurrentExceptionMsg()
        stderr.writeLine "Use --help for usage information"
        quit QuitFailure

proc logArgs*(opts: ref) {.discardable.} =
    if opts.bam.len > 0:
        log("ARG", fmt"Input BAM(s): {opts.bam}")
    else:
        log("ARG", fmt"Input FASTQ1: {opts.fastq1}")
        if opts.fastq2.len > 0:
            log("ARG", fmt"Input FASTQ2: {opts.fastq2}")
        if opts.sample.len > 0:
            log("ARG", fmt"Sample ID: {opts.sample}")
    log("ARG", fmt"Output: {opts.out}")
    log("ARG", fmt"Subset: {opts.subset}")
    log("ARG", fmt"MinQ: {opts.minq}")
    log("ARG", fmt"Write per-read GC: {opts.write_read_gc}")
    log("ARG", fmt"Write run stats: {opts.write_run_stats}")
