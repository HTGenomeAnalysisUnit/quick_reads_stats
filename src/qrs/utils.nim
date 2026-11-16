import times
import strutils
import strformat
import system
import math
import os

proc log*(level: string = "INFO", msg: string) {.discardable.} =
    let t = now()
    let f = initTimeFormat("HH:mm:ss-fff")
    var time_string = format(t, f)
    let log_msg = fmt"[{time_string}] - {level}: {msg}"
    stderr.write(log_msg & "\n")

proc elapsed_time*(start_time: float): string =
    let interval = cpuTime() - start_time
    let s = floor(interval)
    let m = floor(((interval - s) * 1000))
    let time_interval = initDuration(seconds = int(s), milliseconds = int(m))
    result = $time_interval

proc progress_counter*(n: int, interval: var int, t0: var float): (bool, string) {.discardable.} =
    result = (false, "")

    if floorMod(n, interval) == 0:
        result = (true, fmt"{n} reads processed. Last batch: {interval} in {elapsed_time(t0)}")
        t0 = cpuTime()
    case n
        of 1000000: interval = 500000
        of 5000000: interval = 1000000
        of 15000000: interval = 5000000
        else: discard

proc read_list*(list: string, name: string = ""): seq[string] =
    if list == "":
        result = @[]
    else:
        if os.fileExists(list):
            log("INFO", fmt"Reading {name} list from file: {list}")
            for line in lines(list):
                result.add(line)
        else:
            log("INFO", fmt"Reading {name} list from comma-separated string")
            result = list.split(",")

        log("INFO", fmt"Imported {result.len} {name}s")

proc deriveSampleName*(fastqPath: string): string =
    var base = extractFilename(fastqPath)
    let lower = base.toLowerAscii()
    for ext in [".fastq.gz", ".fq.gz", ".fastq.bgz", ".fq.bgz", ".fastq", ".fq"]:
        if lower.endsWith(ext) and base.len > ext.len:
            base.setLen(base.len - ext.len)
            break
    for suffix in ["_R1_001", ".R1_001", "_R1", ".R1", "_1"]:
        if base.endsWith(suffix) and base.len > suffix.len:
            base.setLen(base.len - suffix.len)
            break
    result = base

proc normFastqNameBounds*(raw: cstring, rawLen: int): tuple[s: int, e: int] {.inline.} =
    var endLen = rawLen
    if endLen > 0 and raw[endLen - 1] == '\r':
        dec endLen
    let s = if endLen > 0 and raw[0] == '@': 1 else: 0
    var e = s
    while e < endLen and raw[e] != ' ' and raw[e] != '\t':
        inc e
    if e - s > 2 and raw[e - 2] == '/' and (raw[e - 1] == '1' or raw[e - 1] == '2'):
        e -= 2
    result = (s, e)

proc sameNormQname*(raw1: cstring, s1, e1: int, raw2: cstring, s2, e2: int): bool {.inline.} =
    let len1 = e1 - s1
    if len1 != (e2 - s2):
        return false
    if len1 <= 0:
        return true
    let p1 = cast[pointer](cast[uint](raw1) + uint(s1))
    let p2 = cast[pointer](cast[uint](raw2) + uint(s2))
    result = equalMem(p1, p2, len1)

proc hashNormQname*(raw: cstring, s, e: int): uint64 {.inline.} =
    # 64-bit FNV-1a followed by MurmurHash3 64-bit finalizer mix
    var h: uint64 = 14695981039346656037'u64
    for i in s..<e:
        h = (h xor uint64(ord(raw[i]))) * 1099511628211'u64
    h = h xor (h shr 33)
    h = h * 0xff51afd7ed558ccd'u64
    h = h xor (h shr 33)
    h = h * 0xc4ceb9fe1a85ec53'u64
    h = h xor (h shr 33)
    result = h

proc extractRunKey*(raw: cstring, s, e: int): string {.inline.} =
    # Extract first 4 colon-separated fields (INST:RUN:FLOW:LANE)
    var colons = 0
    var stopIdx = -1
    for i in s..<e:
        if raw[i] == ':':
            inc colons
            if colons == 4:
                stopIdx = i
                break
    if colons == 3:
        stopIdx = e
    if stopIdx <= s or colons < 3:
        return ""
    let keyLen = stopIdx - s
    result = newString(keyLen)
    for i in 0..<keyLen:
        result[i] = raw[s + i]

proc histMean*(h: array[256, int64]): float =
    var total: int64 = 0
    var weighted: int64 = 0
    for i in 0..255:
        let c = h[i]
        if c > 0:
            total += c
            weighted += int64(i) * c
    if total == 0:
        return 0.0
    result = weighted.float / total.float

proc histMedian*(h: array[256, int64]): float =
    var total: int64 = 0
    for i in 0..255:
        total += h[i]
    if total == 0:
        return 0.0
    let mid1 = (total - 1) div 2
    let mid2 = total div 2
    var cum: int64 = 0
    var val1 = -1
    var val2 = -1
    for i in 0..255:
        cum += h[i]
        if val1 < 0 and cum > mid1:
            val1 = i
        if cum > mid2:
            val2 = i
            break
    result = (val1 + val2).float / 2.0

proc histFracAbove*(h: array[256, int64], minVal: int): float =
    var total: int64 = 0
    var above: int64 = 0
    for i in 0..255:
        let c = h[i]
        total += c
        if i >= minVal:
            above += c
    if total == 0:
        return 0.0
    result = above.float / total.float
