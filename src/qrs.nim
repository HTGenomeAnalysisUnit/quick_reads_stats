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
import std/sequtils
import math
import random
from os import createDir, walkPattern, fileExists

const VERSION = "0.2.0"
const READNAMEIDX = (inst: 0, run: 1, flow: 2, lane: 3)

type Region = tuple
    chrom: string
    start: int
    nreads: int

type SampleStats = object
    is_fastq: bool
    is_paired_fastq: bool
    n_reads: int64
    n_aligned_reads: int64
    sum_readlen: int64
    min_readlen: int
    base_q_hist: array[256, int64]
    mq_hist: array[256, int64]
    n_dup: int64
    n_proper_pair: int64
    total_fastq_reads: int64
    paired_fastq_reads: int64
    sum_gc_perc: float64
    sum_error_rate: float64
    sum_indels: int64
    read_names: HashSet[string]
    gc_written_idx: int64

type FastqSampledRead = object
    read_len: int
    gc_perc: float64
    qual: seq[uint8]
    run_key: string

type FastqReader = object
    fp: ptr htsFile
    hdr: kstring_t
    seqStr: kstring_t
    plus: kstring_t
    qual: kstring_t
    eof: bool
    path: string

proc initSampleStats(isFastq: bool = false, isPairedFastq: bool = false): SampleStats =
    result.is_fastq = isFastq
    result.is_paired_fastq = isPairedFastq
    result.min_readlen = high(int)

proc randOpen01(): float {.inline.} =
    result = rand(1.0)
    while result <= 0.0 or result >= 1.0:
        result = rand(1.0)

proc expandInputPaths(inputs: seq[string]): seq[string] =
    for pat in inputs:
        if '*' in pat or '?' in pat or '[' in pat:
            var matched: seq[string]
            for f in walkPattern(pat):
                matched.add(f)
            matched.sort()
            if matched.len > 0:
                result.add(matched)
            else:
                result.add(pat)
        else:
            result.add(pat)

proc computeCigarStats(r: Record): tuple[is_aligned: bool, error_rate: float, indels: int] {.inline.} =
    if r.flag.unmapped:
        return (false, 0.0, 0)
    var
        aligned_qlen = 0
        cigar_diffs = 0
        ins_bases = 0
        del_bases = 0
        indels = 0
    for el in r.cigar:
        let oplen = el.len
        case el.op
        of CigarOp.match, CigarOp.equal:
            aligned_qlen += oplen
        of CigarOp.diff:
            aligned_qlen += oplen
            cigar_diffs += oplen
        of CigarOp.insert:
            aligned_qlen += oplen
            ins_bases += oplen
            inc indels
        of CigarOp.deletion:
            del_bases += oplen
            inc indels
        else:
            discard

    if aligned_qlen <= 0:
        return (false, 0.0, 0)

    let nm_opt = tag[int](r, "NM")
    let edit_dist = if nm_opt.isSome: max(0, nm_opt.get) else: (cigar_diffs + ins_bases + del_bases)
    result = (true, edit_dist.float / aligned_qlen.float, indels)

proc openFastq(r: var FastqReader, path: string) =
    r.path = path
    r.eof = false
    if not fileExists(path):
        log("ERROR", fmt"FASTQ file not found: {path}")
        quit QuitFailure
    r.fp = hts_open(cstring(path), "r")
    if r.fp == nil:
        log("ERROR", fmt"Unable to open FASTQ file: {path}")
        quit QuitFailure

proc closeFastq(r: var FastqReader) =
    if r.hdr.s != nil:
        free(r.hdr.s)
        r.hdr.s = nil
    if r.seqStr.s != nil:
        free(r.seqStr.s)
        r.seqStr.s = nil
    if r.plus.s != nil:
        free(r.plus.s)
        r.plus.s = nil
    if r.qual.s != nil:
        free(r.qual.s)
        r.qual.s = nil
    if r.fp != nil:
        discard hts_close(r.fp)
        r.fp = nil

proc stripCR(k: var kstring_t) {.inline.} =
    if k.l > 0 and k.s[int(k.l) - 1] == '\r':
        dec k.l
        k.s[int(k.l)] = '\0'

proc nextRecord(r: var FastqReader): bool =
    if r.eof:
        return false
    var l1 = hts_getline(r.fp, cint(10), addr r.hdr)
    while l1 == 0:
        l1 = hts_getline(r.fp, cint(10), addr r.hdr)
    if l1 == -1:
        r.eof = true
        return false
    if l1 < -1:
        raise newException(IOError, fmt"Error reading FASTQ stream from {r.path}")

    let l2 = hts_getline(r.fp, cint(10), addr r.seqStr)
    let l3 = hts_getline(r.fp, cint(10), addr r.plus)
    let l4 = hts_getline(r.fp, cint(10), addr r.qual)
    if l2 < 0 or l3 < 0 or l4 < 0:
        raise newException(ValueError, fmt"Truncated FASTQ record in {r.path}")

    stripCR(r.hdr)
    stripCR(r.seqStr)
    stripCR(r.plus)
    stripCR(r.qual)

    if r.hdr.l == 0 or r.hdr.s[0] != '@':
        raise newException(ValueError, fmt"Malformed FASTQ header line in {r.path}: {$r.hdr.s}")
    if r.plus.l == 0 or r.plus.s[0] != '+':
        raise newException(ValueError, fmt"Malformed FASTQ '+' separator line in {r.path}")
    if r.seqStr.l != r.qual.l:
        raise newException(ValueError, fmt"FASTQ sequence and quality length mismatch in {r.path} ({r.seqStr.l} vs {r.qual.l})")

    result = true

proc countGcRaw(seqPtr: cstring, seqLen: int): int {.inline.} =
    for i in 0..<seqLen:
        let ch = seqPtr[i]
        if ch == 'G' or ch == 'C' or ch == 'g' or ch == 'c':
            inc result

proc seqStatsHeader(minQ: int): string =
    &"SAMPLE\tMEAN_READLEN\tMIN_READLEN\tMEAN_BASEQ\tMEDIAN_BASEQ\tPERC_BASES_ABOVE_Q{minQ}\tMEAN_MAPQ\tMEDIAN_MAPQ\tDUP_RATE\tPROPER_PAIR_RATE\tMEAN_READ_GC_PERC\tMEAN_ERROR_RATE\tMEAN_INDELS\n"

proc save_seq_stats(x: OrderedTable[string, SampleStats], minQ: int = 30, outprefix: string) =
    var persample_out = open(fmt"{outprefix}.read_stats.tsv", fmWrite)
    defer: persample_out.close()
    persample_out.write(seqStatsHeader(minQ))

    for sid, sdata in x.pairs:
        if sdata.n_reads == 0:
            continue
        let
            mean_readlen = sdata.sum_readlen.float / sdata.n_reads.float
            min_readlen = if sdata.min_readlen == high(int): 0 else: sdata.min_readlen
            meanq = histMean(sdata.base_q_hist)
            medianq = histMedian(sdata.base_q_hist)
            perc_minq = histFracAbove(sdata.base_q_hist, minQ)
            mean_gc = sdata.sum_gc_perc / sdata.n_reads.float

        var
            mean_mq_str = "NA"
            median_mq_str = "NA"
            dup_rate_str = "NA"
            proper_pair_str = "NA"
            mean_err_str = "NA"
            mean_indels_str = "NA"

        if not sdata.is_fastq:
            mean_mq_str = histMean(sdata.mq_hist).formatFloat(precision = 4)
            median_mq_str = $histMedian(sdata.mq_hist)
            dup_rate_str = (sdata.n_dup.float / sdata.n_reads.float).formatFloat(precision = 4)
            proper_pair_str = (sdata.n_proper_pair.float / sdata.n_reads.float).formatFloat(precision = 4)
            if sdata.n_aligned_reads > 0:
                mean_err_str = (sdata.sum_error_rate / sdata.n_aligned_reads.float).formatFloat(precision = 4)
                mean_indels_str = (sdata.sum_indels.float / sdata.n_aligned_reads.float).formatFloat(precision = 4)
        else:
            if sdata.is_paired_fastq and sdata.total_fastq_reads > 0:
                proper_pair_str = (sdata.paired_fastq_reads.float / sdata.total_fastq_reads.float).formatFloat(precision = 4)

        let result_line = &"{sid}\t{mean_readlen.formatFloat(precision = 4)}\t{min_readlen}\t{meanq.formatFloat(precision = 4)}\t{medianq}\t{perc_minq.formatFloat(precision = 4)}\t{mean_mq_str}\t{median_mq_str}\t{dup_rate_str}\t{proper_pair_str}\t{mean_gc.formatFloat(precision = 4)}\t{mean_err_str}\t{mean_indels_str}\n"
        persample_out.write(result_line)

proc save_run_stats(x: OrderedTable[string, SampleStats], outprefix: string) =
    let header_details = "SAMPLE\tINSTRUMENT\tRUN\tFLOWCELL\tLANE\n"
    let header_persample = "SAMPLE\tN_RUNS\tN_FLOWCELLS\tN_LANES\n"

    var details_out = open(fmt"{outprefix}.rundetails.tsv", fmWrite)
    defer: details_out.close()
    var persample_out = open(fmt"{outprefix}.per_sample_run.tsv", fmWrite)
    defer: persample_out.close()

    details_out.write(header_details)
    persample_out.write(header_persample)

    for sid, sdata in x.pairs:
        if sdata.n_reads == 0:
            continue
        var
            run_set: HashSet[string]
            flow_set: HashSet[string]
            n_lanes = 0
        for read_name in sdata.read_names:
            let toks = read_name.split(":")
            if toks.len < 4:
                continue
            inc n_lanes
            run_set.incl(fmt"{toks[READNAMEIDX.inst]}:{toks[READNAMEIDX.run]}")
            flow_set.incl(fmt"{toks[READNAMEIDX.inst]}:{toks[READNAMEIDX.run]}:{toks[READNAMEIDX.flow]}")
            details_out.write(&"{sid}\t{toks[READNAMEIDX.inst]}\t{toks[READNAMEIDX.run]}\t{toks[READNAMEIDX.flow]}\t{toks[READNAMEIDX.lane]}\n")

        persample_out.write(&"{sid}\t{run_set.len}\t{flow_set.len}\t{n_lanes}\n")

proc processBamInputs(bamInputs: seq[string], max_reads: int, write_gc: bool, write_run: bool, outprefix: string, samplemeta: var OrderedTable[string, SampleStats]) =
    let bamPaths = expandInputPaths(bamInputs)
    let gcOutdir = fmt"{outprefix}.GC_stats"
    if write_gc:
        createDir(gcOutdir)

    var
        b: Bam
        b_quals: seq[uint8]
        read_seq: string
        gcInitialized: HashSet[string]

    let reg_exp = re"(chr){0,1}[0-9X]+$"

    for inbam in bamPaths:
        echo fmt"Start reading {inbam}"
        if not open(b, cstring(inbam), index = true):
            log("ERROR", fmt"Failed to open BAM/CRAM (or its index): {inbam}")
            quit QuitFailure

        var
            n_mapped = 0
            tot_length = 0
            random_idx: seq[Region]
            rnd_region: Region

        let targets = b.hdr.targets.filterIt(it.name.match(reg_exp))
        for t in targets:
            n_mapped += b.idx.stats(t.tid).mapped.int
            tot_length += t.length.int

        if max_reads == -1 or n_mapped < max_reads:
            for t in targets:
                rnd_region.chrom = t.name
                rnd_region.start = 1
                rnd_region.nreads = b.idx.stats(t.tid).mapped.int
                random_idx.add(rnd_region)
        else:
            for t in targets:
                let target_weight = t.length.int / tot_length
                let rnd_per_target = int(max_reads.float * target_weight)
                let max_start = max(0, b.idx.stats(t.tid).mapped.int - rnd_per_target)

                rnd_region.chrom = t.name
                rnd_region.start = rand(max_start)
                rnd_region.nreads = rnd_per_target
                random_idx.add(rnd_region)

        var rg_def: Table[string, string]
        var sample_ids: HashSet[string]
        var mysample = ""
        let header_str = $b.hdr
        for line in header_str.split("\n"):
            if line.startsWith("@RG"):
                var tag_tb: Table[string, string]
                let tags = line.split("\t")
                for t in tags[1..tags.high]:
                    let x = t.split(":")
                    if x.len >= 2:
                        tag_tb[x[0]] = x[1]
                if "ID" in tag_tb and "SM" in tag_tb:
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
            b.close()
            continue

        if mysample == "":
            mysample = deriveSampleName(inbam)

        if mysample notin samplemeta:
            samplemeta[mysample] = initSampleStats(isFastq = false)

        var gcFile: File
        var gcOpen = false
        if write_gc:
            let gcPath = fmt"{gcOutdir}/{mysample}.GC_stats.tsv"
            if mysample notin gcInitialized:
                gcFile = open(gcPath, fmWrite)
                gcFile.write("READ_IDX\tGC_PERC\n")
                gcInitialized.incl(mysample)
            else:
                gcFile = open(gcPath, fmAppend)
            gcOpen = true

        if n_mapped == 0:
            echo fmt"WARNING: No mapped reads found in {inbam}"
            if gcOpen:
                gcFile.close()
            b.close()
            continue

        var
            t0 = cpuTime()
            n = 0
            interval = 100000

        template sdata(): var SampleStats = samplemeta[mysample]

        for region in random_idx:
            if region.nreads <= 0:
                continue
            var n_per_target = 0
            for r in b.query(region.chrom, region.start):
                if n_per_target == region.nreads:
                    break
                inc n_per_target
                inc n
                var (dolog, log_msg) = progress_counter(n, interval, t0)
                if dolog:
                    log("INFO", log_msg)

                read_seq = r.sequence(read_seq)
                b_quals = r.base_qualities(b_quals)
                let
                    read_len = read_seq.len
                    gc_n = read_seq.count({'G', 'C', 'g', 'c'})
                    gc_perc = if read_len > 0: gc_n.float / read_len.float else: 0.0

                inc sdata.n_reads
                sdata.sum_readlen += read_len
                if read_len < sdata.min_readlen:
                    sdata.min_readlen = read_len
                sdata.mq_hist[int(r.mapping_quality)] += 1
                for q in b_quals:
                    sdata.base_q_hist[int(q)] += 1
                if r.flag.dup:
                    inc sdata.n_dup
                if r.flag.proper_pair:
                    inc sdata.n_proper_pair
                sdata.sum_gc_perc += gc_perc

                let cigStats = computeCigarStats(r)
                if cigStats.is_aligned:
                    inc sdata.n_aligned_reads
                    sdata.sum_error_rate += cigStats.error_rate
                    sdata.sum_indels += cigStats.indels

                if write_run:
                    let qn = r.qname
                    let (qs, qe) = normFastqNameBounds(cstring(qn), qn.len)
                    let rkey = extractRunKey(cstring(qn), qs, qe)
                    if rkey.len > 0:
                        sdata.read_names.incl(rkey)

                if gcOpen:
                    gcFile.write(&"{sdata.gc_written_idx}\t{gc_perc.formatFloat(precision = 4)}\n")
                    inc sdata.gc_written_idx

        if gcOpen:
            gcFile.close()
        b.close()

proc fillSampledRead(dst: var FastqSampledRead, r: FastqReader, write_run: bool) {.inline.} =
    let seqLen = int(r.seqStr.l)
    let gc_n = countGcRaw(r.seqStr.s, seqLen)
    dst.read_len = seqLen
    dst.gc_perc = if seqLen > 0: gc_n.float / seqLen.float else: 0.0
    dst.qual.setLen(seqLen)
    for i in 0..<seqLen:
        let q = clamp(ord(r.qual.s[i]) - 33, 0, 255)
        dst.qual[i] = uint8(q)
    if write_run:
        let (qs, qe) = normFastqNameBounds(r.hdr.s, int(r.hdr.l))
        dst.run_key = extractRunKey(r.hdr.s, qs, qe)

proc streamFastqRecordDirect(sdata: var SampleStats, r: FastqReader, write_run: bool, gcOpen: bool, gcFile: var File) {.inline.} =
    let seqLen = int(r.seqStr.l)
    let gc_n = countGcRaw(r.seqStr.s, seqLen)
    let gc_perc = if seqLen > 0: gc_n.float / seqLen.float else: 0.0

    inc sdata.n_reads
    sdata.sum_readlen += seqLen
    if seqLen < sdata.min_readlen:
        sdata.min_readlen = seqLen
    for i in 0..<seqLen:
        let q = clamp(ord(r.qual.s[i]) - 33, 0, 255)
        sdata.base_q_hist[q] += 1
    sdata.sum_gc_perc += gc_perc

    if write_run:
        let (qs, qe) = normFastqNameBounds(r.hdr.s, int(r.hdr.l))
        let rkey = extractRunKey(r.hdr.s, qs, qe)
        if rkey.len > 0:
            sdata.read_names.incl(rkey)

    if gcOpen:
        gcFile.write(&"{sdata.gc_written_idx}\t{gc_perc.formatFloat(precision = 4)}\n")
        inc sdata.gc_written_idx

proc processFastqInputs(fq1Path: string, fq2Path: string, sampleArg: string, max_reads: int, write_gc: bool, write_run: bool, outprefix: string, samplemeta: var OrderedTable[string, SampleStats]) =
    let isPaired = fq2Path.len > 0
    let mysample = if sampleArg.len > 0: sampleArg else: deriveSampleName(fq1Path)

    echo "=== FASTQ INFO ==="
    echo fmt"Sample: {mysample}"
    echo fmt"FASTQ1: {fq1Path}"
    if isPaired:
        echo fmt"FASTQ2: {fq2Path}"
    echo "======================"

    var sdata = initSampleStats(isFastq = true, isPairedFastq = isPaired)

    var fq1: FastqReader
    var fq2: FastqReader
    openFastq(fq1, fq1Path)
    defer: closeFastq(fq1)
    if isPaired:
        openFastq(fq2, fq2Path)
    defer:
        if isPaired:
            closeFastq(fq2)

    let gcOutdir = fmt"{outprefix}.GC_stats"
    var gcFile: File
    var gcOpen = false
    if write_gc:
        createDir(gcOutdir)
        gcFile = open(fmt"{gcOutdir}/{mysample}.GC_stats.tsv", fmWrite)
        gcFile.write("READ_IDX\tGC_PERC\n")
        gcOpen = true
    defer:
        if gcOpen:
            gcFile.close()

    var
        reservoir: seq[FastqSampledRead]
        streamIdx: int64 = 0
        algL_w: float = 0.0
        nextSampleIdx: int64 = 0
        pendingHashes: HashSet[uint64]
        t0 = cpuTime()
        n = 0
        interval = 100000

    if max_reads > 0:
        reservoir = newSeqOfCap[FastqSampledRead](max_reads)

    template handleRecord(reader: FastqReader) =
        inc n
        var (dolog, log_msg) = progress_counter(n, interval, t0)
        if dolog:
            log("INFO", log_msg)

        if max_reads == -1:
            streamFastqRecordDirect(sdata, reader, write_run, gcOpen, gcFile)
        else:
            if streamIdx < int64(max_reads):
                var item: FastqSampledRead
                fillSampledRead(item, reader, write_run)
                reservoir.add(item)
                if streamIdx == int64(max_reads - 1):
                    algL_w = exp(ln(randOpen01()) / max_reads.float)
                    let denom = ln(1.0 - algL_w)
                    let skip = if denom < 0.0: int64(floor(ln(randOpen01()) / denom)) else: 0'i64
                    nextSampleIdx = int64(max_reads) + skip
            elif streamIdx == nextSampleIdx:
                let slot = rand(max_reads - 1)
                fillSampledRead(reservoir[slot], reader, write_run)
                algL_w = algL_w * exp(ln(randOpen01()) / max_reads.float)
                let denom = ln(1.0 - algL_w)
                let skip = if denom < 0.0: int64(floor(ln(randOpen01()) / denom)) else: 0'i64
                nextSampleIdx = streamIdx + 1'i64 + skip
            inc streamIdx

    if not isPaired:
        while fq1.nextRecord():
            inc sdata.total_fastq_reads
            handleRecord(fq1)
    else:
        while true:
            let got1 = fq1.nextRecord()
            let got2 = fq2.nextRecord()
            if not got1 and not got2:
                break

            if got1 and got2:
                sdata.total_fastq_reads += 2
                let (s1, e1) = normFastqNameBounds(fq1.hdr.s, int(fq1.hdr.l))
                let (s2, e2) = normFastqNameBounds(fq2.hdr.s, int(fq2.hdr.l))
                if sameNormQname(fq1.hdr.s, s1, e1, fq2.hdr.s, s2, e2):
                    sdata.paired_fastq_reads += 2
                else:
                    let h1 = hashNormQname(fq1.hdr.s, s1, e1)
                    if h1 in pendingHashes:
                        pendingHashes.excl(h1)
                        sdata.paired_fastq_reads += 2
                    else:
                        pendingHashes.incl(h1)

                    let h2 = hashNormQname(fq2.hdr.s, s2, e2)
                    if h2 in pendingHashes:
                        pendingHashes.excl(h2)
                        sdata.paired_fastq_reads += 2
                    else:
                        pendingHashes.incl(h2)

                handleRecord(fq1)
                handleRecord(fq2)
            elif got1:
                inc sdata.total_fastq_reads
                if pendingHashes.len > 0:
                    let (s1, e1) = normFastqNameBounds(fq1.hdr.s, int(fq1.hdr.l))
                    let h1 = hashNormQname(fq1.hdr.s, s1, e1)
                    if h1 in pendingHashes:
                        pendingHashes.excl(h1)
                        sdata.paired_fastq_reads += 2
                handleRecord(fq1)
            else:
                inc sdata.total_fastq_reads
                if pendingHashes.len > 0:
                    let (s2, e2) = normFastqNameBounds(fq2.hdr.s, int(fq2.hdr.l))
                    let h2 = hashNormQname(fq2.hdr.s, s2, e2)
                    if h2 in pendingHashes:
                        pendingHashes.excl(h2)
                        sdata.paired_fastq_reads += 2
                handleRecord(fq2)

    if max_reads > 0 and reservoir.len > 0:
        for item in reservoir:
            inc sdata.n_reads
            sdata.sum_readlen += item.read_len
            if item.read_len < sdata.min_readlen:
                sdata.min_readlen = item.read_len
            for q in item.qual:
                sdata.base_q_hist[int(q)] += 1
            sdata.sum_gc_perc += item.gc_perc
            if write_run and item.run_key.len > 0:
                sdata.read_names.incl(item.run_key)
            if gcOpen:
                gcFile.write(&"{sdata.gc_written_idx}\t{item.gc_perc.formatFloat(precision = 4)}\n")
                inc sdata.gc_written_idx

    if sdata.n_reads == 0:
        echo fmt"WARNING: No reads found in FASTQ input for {mysample}"
    else:
        echo fmt"Total FASTQ reads scanned: {sdata.total_fastq_reads} (sampled: {sdata.n_reads})"

    samplemeta[mysample] = sdata

proc main*() =
    echo fmt"Start Quick Reads Stats - VERSION {VERSION}"
    randomize()
    var opts = parseCmdLine()
    opts.logArgs()

    let max_reads = parseInt(opts.subset)
    let minq = parseInt(opts.minq)
    let start_time = cpuTime()
    var samplemeta: OrderedTable[string, SampleStats]

    if opts.bam.len > 0:
        processBamInputs(opts.bam, max_reads, opts.write_read_gc, opts.write_run_stats, opts.out, samplemeta)
    else:
        processFastqInputs(opts.fastq1, opts.fastq2, opts.sample, max_reads, opts.write_read_gc, opts.write_run_stats, opts.out, samplemeta)

    echo "=== SAVING OUTPUT STATS ==="
    var t0 = cpuTime()
    if opts.write_run_stats:
        samplemeta.save_run_stats(opts.out)
        echo fmt"Run stats saved in {elapsed_time(t0)}"
        t0 = cpuTime()

    samplemeta.save_seq_stats(minq, opts.out)
    echo fmt"Sequencing stats saved in {elapsed_time(t0)}"

    if opts.write_read_gc:
        echo "Per-read GC stats streamed to disk"

    echo fmt"Completed in {elapsed_time(start_time)}"

when isMainModule:
    main()
