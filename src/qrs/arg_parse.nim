import argparse
import strformat
from ./utils import log

var p = newParser("CSQ Selector"):
    help("Compute Illumina run stats from read names")
    arg("bam", nargs = -1, help="input BAM/CRAM file(s). Glob pattern allowed")
    option("-o", "--out", help="Output file", default=some("run_stats.tsv"))
    option("-s", "--subset", help="Read only N reads from each file. Use -1 for all reads", default=some("1000000"))
    option("-q", "--minq", help="Min base quality for stats reporting", default=some("30"))
        
proc parseCmdLine*(): ref =
    try:
        result = p.parse() 
    except ShortCircuit as e:
        if e.flag == "argparse_help":
            echo p.help
            quit QuitSuccess
    except UsageError:
        stderr.writeLine getCurrentExceptionMsg() 
        echo "Use --help for usage information"
        quit QuitSuccess

proc logArgs*(opts: ref) {.discardable.} =
    log("ARG", fmt"Input file: {opts.bam}")
    log("ARG", fmt"Output: {opts.out}")
    log("ARG", fmt"Subset: {opts.subset}")