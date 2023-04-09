# This is just an example to get you started. A typical binary package
# uses this file as the main entry point of the application.

import hts
import re
import strutils
import qrs/arg_parse
import qrs/utils
import sets
import tables
import times
import strformat
import std/algorithm
import stats
import math
import random
from os import createDir

const VERSION="0.1.1"
const READNAMEIDX = (inst: 0, run: 1, flow: 2, lane: 3)

type Region = tuple
    chrom: string
    start: int
    nreads: int

type ReadData = object
    read_names: HashSet[string]
    read_len: seq[int]
    mq: seq[uint8]
    gc_perc: seq[float]
    base_q: seq[uint8]
    dup: seq[int]

# Calculate the median of a sequence of integers
proc median(numbers: seq[int]): float =
  # Sort the sequence
  var s = numbers
  s.sort()

  # Get the length of the sequence
  let n = s.len

  # If the length of the sequence is even, return the average of the
  # two middle elements
  if floorMod(n, 2) == 0:
    return (s[n div 2 - 1] + s[n div 2]) / 2
  # If the length of the sequence is odd, return the middle element
  else:
    return s[n div 2].float

    
proc toInt(x: seq[uint8]): seq[int] {.inline.} =
    for val in x: result.add(int(val))

proc count_above(x:seq[int], min: int): int =
    result = 0
    for val in x:
        if val >= min: result += 1

proc save_gc_stats(x: Table[string, ReadData], outprefix: string) =
    let header_persample = &"READ_IDX\tGC_PERC\n"
    let outdir = fmt"{outprefix}.GC_stats"
    createDir(outdir)
    for sid, sdata in x.pairs():
        var persample_out = open(fmt"{outdir}/{sid}.GC_stats.tsv", fmWrite)
        persample_out.write(header_persample)
        for i in 0..sdata.gc_perc.high:
            let result_line = &"{i}\t{sdata.gc_perc[i].formatFloat(precision = 4)}\n"
            persample_out.write(result_line)
        persample_out.close()

proc save_seq_stats(x: Table[string, ReadData], minQ: int = 30, outprefix: string) =
    let header_persample = &"SAMPLE\tMEAN_READLEN\tMIN_READLEN\tMEAN_BASEQ\tMEDIAN_BASEQ\tPERC_BASES_ABOVE_Q{minQ}\tMEAN_MAPQ\tMEDIAN_MAPQ\tDUP_RATE\tMEAN_READ_GC_PERC\n"
    
    var persample_out = open(fmt"{outprefix}.reads_stat.tsv", fmWrite)
    persample_out.write(header_persample)
    for sid, sdata in x.pairs:
        let 
            base_q_int = sdata.base_q.toInt
            map_q_int = sdata.mq.toInt
            meanq = mean(base_q_int)
            medianq = median(base_q_int)
            perc_minq = base_q_int.count_above(minQ) / base_q_int.len
            mean_mq = mean(map_q_int)
            median_mq = median(map_q_int)
            dup_rate = sum(sdata.dup) / sdata.dup.len
            mean_gc = mean(sdata.gc_perc)
            mean_readlen = mean(sdata.read_len)
            min_readlen = min(sdata.read_len)
        
        let result_line = &"{sid}\t{mean_readlen.formatFloat(precision = 4)}\t{min_readlen}\t{meanq.formatFloat(precision = 4)}\t{medianq}\t{perc_minq.formatFloat(precision = 4)}\t{mean_mq.formatFloat(precision = 4)}\t{median_mq}\t{dup_rate.formatFloat(precision = 4)}\t{mean_gc.formatFloat(precision = 4)}\n"
        persample_out.write(result_line)
    persample_out.close()

proc save_run_stats(x: Table[string, ReadData], outprefix: string) =
    let header_details = "SAMPLE\tINSTRUMENT\tRUN\tFLOWCELL\tLANE\n"
    let header_persample = "SAMPLE\tN_RUNS\tN_FLOWCELLS\tN_LANES\n"
    
    var details_out = open(fmt"{outprefix}.rundetails.tsv", fmWrite)
    var persample_out = open(fmt"{outprefix}.per_sample_run.tsv", fmWrite)
    details_out.write(header_details)
    persample_out.write(header_persample)
    for sid, sdata in x.pairs:
        var 
            run_set: HashSet[string]
            flow_set: HashSet[string]
            n_lanes: int
        for read_name in sdata.read_names:
            let toks = read_name.split(":")
            n_lanes += 1
            run_set.incl(fmt"{toks[READNAMEIDX.inst]}:{toks[READNAMEIDX.run]}")
            flow_set.incl(fmt"{toks[READNAMEIDX.inst]}:{toks[READNAMEIDX.run]}:{toks[READNAMEIDX.flow]}")
            details_out.write(&"{sid}\t{toks[READNAMEIDX.inst]}\t{toks[READNAMEIDX.run]}\t{toks[READNAMEIDX.flow]}\t{toks[READNAMEIDX.lane]}\n")
        
        persample_out.write(&"{sid}\t{run_set.len}\t{flow_set.len}\t{n_lanes}\n")
    
    details_out.close()
    persample_out.close()

proc main*() =
    echo fmt"Start Quick Reads Stats - VERSION {VERSION}"
    var opts = parseCmdLine()
    opts.logArgs()
    
    let max_reads = parseInt(opts.subset)
    let minq = parseInt(opts.minq)
    var 
        samplemeta: Table[string,ReadData]
        b:Bam
        start_time = cpuTime()
        emptyset: ReadData
        b_quals: seq[uint8]
        read_seq: string
        mysample:string

    for inbam in opts.bam:
        echo fmt"Start reading {inbam}"
        open(b, cstring(inbam), index=true)

        #Set random N reads to process on canonical chrom
        var 
            n_mapped = 0
            tot_length = 0
            random_idx: seq[Region]
            rnd_region: Region
            targets = b.hdr.targets
            t_idx = 0
        let reg_exp = re"(chr){0,1}[0-9X]+$"
        for t in targets:
            if not t.name.match(reg_exp): targets.delete(t_idx)
            t_idx += 1
        for t in targets: 
            n_mapped += b.idx.stats(t.tid).mapped.int
            tot_length += t.length.int

        if max_reads == -1 or n_mapped < max_reads:
            #Set the whole target as region to process all reads
            for t in targets:
                rnd_region.chrom = t.name
                rnd_region.start = 1
                rnd_region.nreads = b.idx.stats(t.tid).mapped.int
                random_idx.add(rnd_region)
        else:
            #Get proportion of read to allocate per target
            
            for t in targets:
                let target_weight = t.length.int / tot_length 
                let rnd_per_target = int(max_reads.float * target_weight)

                rnd_region.chrom = t.name
                rnd_region.start = rand((b.idx.stats(t.tid).mapped.int - rnd_per_target))
                rnd_region.nreads = rnd_per_target
                random_idx.add(rnd_region)

        var rg_def: Table[string, string]
        var sample_ids: HashSet[string]
        #Parse sample from header
        let header_str = $b.hdr
        for line in header_str.split("\n"):
            if line.startsWith("@RG"):
                var tag_tb: Table[string, string]
                let tags = line.split("\t")
                for t in tags[1..tags.high]: 
                    let x = t.split(":")
                    tag_tb[x[0]] = x[1]
                rg_def[tag_tb["ID"]] = tag_tb["SM"]
        
        echo "=== BAM INFO ==="
        echo fmt"Mapped reads: {n_mapped}"
        for key, val in rg_def.pairs():
            echo fmt"RG: {key}, Sample: {val}"
            sample_ids.incl(val)
            mysample = val
        echo "======================"
        if sample_ids.len > 1:
            log("WARN", fmt"More than one sample detected in {inbam}. Skipping this file")
            continue
        samplemeta[mysample] = emptyset

        var
            t0 = cpuTime()
            n = 0
            interval = 100000
        
        for region in random_idx:
            var n_per_target = 0
            for r in b.query(region.chrom, region.start):
                if n_per_target == region.nreads: break
                n_per_target += 1
                n += 1
                var (dolog, log_msg) = progress_counter(n, interval, t0)
                if dolog: log("INFO", log_msg)
                
                read_seq = r.sequence(read_seq)
                b_quals = r.base_qualities(b_quals)
                let 
                    toks = r.qname.split(':')
                    read_simple_name = fmt"{toks[0]}:{toks[1]}:{toks[2]}:{toks[3]}"
                    is_dup = (if r.flag.dup: 1 else: 0)
                    read_len = read_seq.len
                    gc_n = read_seq.count({'G','C'}) 

                samplemeta[mysample].read_names.incl(read_simple_name)
                samplemeta[mysample].read_len.add(read_len)
                samplemeta[mysample].mq.add(r.mapping_quality)
                samplemeta[mysample].gc_perc.add(gc_n/read_len)
                samplemeta[mysample].base_q.add(b_quals)
                samplemeta[mysample].dup.add(is_dup)
            
        b.close()
    
    echo "=== MAKING DATAFRAME ==="
    var t0 = cpuTime()
    samplemeta.save_run_stats(opts.out)
    echo fmt"Run stats saved in {elapsed_time(t0)}"
    samplemeta.save_seq_stats(minq, opts.out)
    echo fmt"Sequencing stats saved in {elapsed_time(t0)}"
    samplemeta.save_gc_stats(opts.out)
    echo fmt"GC stats saved in {elapsed_time(t0)}"
    #samplemeta.plotGC(opts.out)
    #echo fmt"GC plots saved in {elapsed_time(t0)}"

    echo fmt"Completed in {elapsed_time(start_time)}"

when isMainModule:
    main()
