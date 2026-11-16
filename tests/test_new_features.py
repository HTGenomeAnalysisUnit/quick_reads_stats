"""Tests for v0.2.0 features: CIGAR error/indel metrics, FASTQ single/paired, GC streaming, and CLI flags."""

import csv
import gzip
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


def qrs_bin():
    return os.environ.get("QRS_BIN", "qrs")


def samtools_bin():
    return os.environ.get("SAMTOOLS", "samtools")


class NewFeaturesTest(unittest.TestCase):
    def test_bam_cigar_error_indels_and_proper_pairs(self):
        with tempfile.TemporaryDirectory(prefix="qrs-cigar-") as tmp:
            root = Path(tmp)
            sam = root / "cigar.sam"
            bam = root / "cigar.bam"
            prefix = root / "out"

            header = [
                "@HD\tVN:1.6\tSO:coordinate",
                "@SQ\tSN:chr1\tLN:10000",
                "@RG\tID:rg1\tSM:cigar_sample",
            ]
            # 4 reads on chr1:
            # 1) 10M, NM:i:1, flag 99 (proper_pair=True, dup=False): aligned=10, err=0.1, indels=0
            # 2) 4M2I4M, NM:i:3, flag 163 (proper_pair=True, dup=False): aligned=10, err=0.3, indels=1
            # 3) 5M3D5M2S (seq len 12), NM:i:4, flag 65 (proper_pair=False, dup=False): aligned=10, err=0.4, indels=1
            # 4) 4=2X1I1D3= (seq len 10, no NM tag), flag 1123 (proper_pair=True, dup=True): aligned=10, err=0.4, indels=2
            records = [
                "A01:10:FC1:1:1:1:1\t99\tchr1\t10\t60\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1\tNM:i:1",
                "A01:10:FC1:1:1:1:2\t163\tchr1\t30\t50\t4M2I4M\t=\t130\t100\tACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1\tNM:i:3",
                "A01:10:FC1:2:1:1:3\t65\tchr1\t60\t40\t5M3D5M2S\t=\t160\t100\tACGTACGTACGT\tIIIIIIIIIIII\tRG:Z:rg1\tNM:i:4",
                "A01:11:FC2:1:1:1:4\t1123\tchr1\t90\t30\t4=2X1I1D3=\t=\t190\t100\tGGGGGGGGGG\t5555555555\tRG:Z:rg1",
            ]
            sam.write_text("\n".join(header + records) + "\n")
            subprocess.run([samtools_bin(), "view", "-b", "-o", str(bam), str(sam)],
                           check=True, capture_output=True)
            subprocess.run([samtools_bin(), "index", str(bam)],
                           check=True, capture_output=True)

            run = subprocess.run(
                [qrs_bin(), "--bam", str(bam), "--subset", "-1",
                 "--write_read_gc", "--write_run_stats", "--out", str(prefix)],
                text=True, capture_output=True,
            )
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)

            with (root / "out.read_stats.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["SAMPLE"], "cigar_sample")
            self.assertAlmostEqual(float(row["MEAN_READLEN"]), (10 + 10 + 12 + 10) / 4, places=3)
            self.assertEqual(int(row["MIN_READLEN"]), 10)
            self.assertAlmostEqual(float(row["MEAN_MAPQ"]), 45.0, places=3)
            self.assertAlmostEqual(float(row["MEDIAN_MAPQ"]), 45.0, places=3)
            self.assertAlmostEqual(float(row["DUP_RATE"]), 0.25, places=4)
            self.assertAlmostEqual(float(row["PROPER_PAIR_RATE"]), 0.75, places=4)
            self.assertAlmostEqual(float(row["MEAN_ERROR_RATE"]), 0.3, places=4)
            self.assertAlmostEqual(float(row["MEAN_INDELS"]), 1.0, places=4)

            with (root / "out.per_sample_run.tsv").open() as handle:
                run_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(run_rows), 1)
            self.assertEqual(int(run_rows[0]["N_RUNS"]), 2)
            self.assertEqual(int(run_rows[0]["N_FLOWCELLS"]), 2)
            self.assertEqual(int(run_rows[0]["N_LANES"]), 3)

    def test_single_end_fastq_gz(self):
        with tempfile.TemporaryDirectory(prefix="qrs-fq-se-") as tmp:
            root = Path(tmp)
            fq1 = root / "mysample_R1_001.fastq.gz"
            prefix = root / "se_stats"

            # 3 reads:
            # Read 0: len 4, GC=2/4=0.5, quals 'IIII' (40,40,40,40)
            # Read 1: len 6, GC=6/6=1.0 (lowercase gc mix), quals '555555' (20,20,20,20,20,20)
            # Read 2: len 4, GC=0/4=0.0, quals 'IIII' (40,40,40,40)
            lines = [
                "@INST1:1:FC1:1:1101:10:20 1:N:0:ATCG",
                "ACGT",
                "+",
                "IIII",
                "@INST1:1:FC1:1:1101:10:21 1:N:0:ATCG",
                "GGccGG",
                "+",
                "555555",
                "@INST1:1:FC1:2:1101:10:22 1:N:0:ATCG",
                "AATT",
                "+",
                "IIII",
            ]
            with gzip.open(fq1, "wt") as gz:
                gz.write("\n".join(lines) + "\n")

            run = subprocess.run(
                [qrs_bin(), "--fastq1", str(fq1), "--subset", "-1",
                 "--write_read_gc", "--out", str(prefix)],
                text=True, capture_output=True,
            )
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertFalse((root / "se_stats.rundetails.tsv").exists())

            with (root / "se_stats.read_stats.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["SAMPLE"], "mysample")
            self.assertAlmostEqual(float(row["MEAN_READLEN"]), 14 / 3, places=3)
            self.assertEqual(int(row["MIN_READLEN"]), 4)
            # 8 bases at Q40, 6 bases at Q20 -> mean = (320 + 120)/14 = 31.43, median = 40.0, frac>=30 = 8/14 = 0.5714
            self.assertAlmostEqual(float(row["MEAN_BASEQ"]), 440 / 14, places=2)
            self.assertAlmostEqual(float(row["MEDIAN_BASEQ"]), 40.0, places=3)
            self.assertAlmostEqual(float(row["PERC_BASES_ABOVE_Q30"]), 8 / 14, places=3)
            self.assertAlmostEqual(float(row["MEAN_READ_GC_PERC"]), 0.5, places=4)
            for na_col in ("MEAN_MAPQ", "MEDIAN_MAPQ", "DUP_RATE", "PROPER_PAIR_RATE", "MEAN_ERROR_RATE", "MEAN_INDELS"):
                self.assertEqual(row[na_col], "NA", f"Expected NA for {na_col}")

            with (root / "se_stats.GC_stats/mysample.GC_stats.tsv").open() as handle:
                gc_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([int(r["READ_IDX"]) for r in gc_rows], [0, 1, 2])
            self.assertEqual([float(r["GC_PERC"]) for r in gc_rows], [0.5, 1.0, 0.0])

    def test_paired_fastq_proper_pairs_resync_and_subsampling(self):
        with tempfile.TemporaryDirectory(prefix="qrs-fq-pe-") as tmp:
            root = Path(tmp)
            fq1 = root / "pe_R1.fastq.gz"
            fq2 = root / "pe_R2.fastq.gz"

            # Create 200 reads in R1 (0..199) and 196 reads in R2 (dropping indices 5, 50, 150, 199)
            # Total reads = 396; paired reads = 196 * 2 = 392 -> PROPER_PAIR_RATE = 392 / 396 = 0.9899
            dropped = {5, 50, 150, 199}
            with gzip.open(fq1, "wt") as g1, gzip.open(fq2, "wt") as g2:
                for i in range(200):
                    lane = 1 if i < 100 else 2
                    g1.write(f"@A01245:148:HMT2:L{lane}:1101:{i}:1000/1\nACGTACGT\n+\nIIIIIIII\n")
                    if i not in dropped:
                        g2.write(f"@A01245:148:HMT2:L{lane}:1101:{i}:1000/2\nGGCCAAtt\n+\nIIIIIIII\n")

            prefix = root / "pe_sub"
            run = subprocess.run(
                [qrs_bin(), "--fastq1", str(fq1), "--fastq2", str(fq2),
                 "--sample", "paired_sample", "--subset", "80",
                 "--write_read_gc", "--write_run_stats", "--out", str(prefix)],
                text=True, capture_output=True,
            )
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)

            with (root / "pe_sub.read_stats.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            row = rows[0]
            self.assertEqual(row["SAMPLE"], "paired_sample")
            self.assertAlmostEqual(float(row["PROPER_PAIR_RATE"]), 392 / 396, places=3)
            self.assertEqual(row["MEAN_ERROR_RATE"], "NA")
            self.assertEqual(row["MEAN_INDELS"], "NA")

            with (root / "pe_sub.GC_stats/paired_sample.GC_stats.tsv").open() as handle:
                gc_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(gc_rows), 80)
            self.assertEqual([int(r["READ_IDX"]) for r in gc_rows], list(range(80)))

            with (root / "pe_sub.per_sample_run.tsv").open() as handle:
                run_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(run_rows), 1)
            self.assertEqual(int(run_rows[0]["N_RUNS"]), 1)
            self.assertEqual(int(run_rows[0]["N_FLOWCELLS"]), 1)
            self.assertEqual(int(run_rows[0]["N_LANES"]), 2)

    def test_cli_validation_errors(self):
        # No input
        r1 = subprocess.run([qrs_bin()], text=True, capture_output=True)
        self.assertNotEqual(r1.returncode, 0)
        # Both --bam and --fastq1
        r2 = subprocess.run([qrs_bin(), "--bam", "a.bam", "--fastq1", "a.fq"],
                            text=True, capture_output=True)
        self.assertNotEqual(r2.returncode, 0)
        # --fastq2 without --fastq1
        r3 = subprocess.run([qrs_bin(), "--fastq2", "b.fq"],
                            text=True, capture_output=True)
        self.assertNotEqual(r3.returncode, 0)

    def test_bam_multiwindow_chromosome_wide_spread(self):
        with tempfile.TemporaryDirectory(prefix="qrs-bam-spread-") as tmp:
            root = Path(tmp)
            sam = root / "spread.sam"
            bam = root / "spread.bam"
            prefix = root / "spread_out"

            # 10 Mb chromosome with 1,000 mapped reads spaced every 10 kb:
            # First 500 reads (0..5 Mb) have GC = 0.0; second 500 reads (5..10 Mb) have GC = 1.0.
            header = [
                "@HD\tVN:1.6\tSO:coordinate",
                "@SQ\tSN:chr1\tLN:10000000",
                "@RG\tID:rg1\tSM:spread_sample",
            ]
            records = []
            for i in range(1000):
                pos = i * 10000 + 100
                seq = "ATATATATAT" if i < 500 else "GCGCGCGCGC"
                records.append(
                    f"A01:1:FC1:L{i}:1:1:1\t0\tchr1\t{pos}\t60\t10M\t*\t0\t0\t{seq}\tIIIIIIIIII\tRG:Z:rg1\tNM:i:0"
                )
            sam.write_text("\n".join(header + records) + "\n")
            subprocess.run([samtools_bin(), "view", "-b", "-o", str(bam), str(sam)],
                           check=True, capture_output=True)
            subprocess.run([samtools_bin(), "index", str(bam)],
                           check=True, capture_output=True)

            run = subprocess.run(
                [qrs_bin(), "--bam", str(bam), "--subset", "200",
                 "--write_read_gc", "--write_run_stats", "--out", str(prefix)],
                text=True, capture_output=True,
            )
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)

            with (root / "spread_out.GC_stats/spread_sample.GC_stats.tsv").open() as handle:
                gc_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(gc_rows), 200)

            with (root / "spread_out.per_sample_run.tsv").open() as handle:
                run_rows = list(csv.DictReader(handle, delimiter="\t"))
            # Every read had a distinct lane tag L{i}, so N_LANES == 200 proves zero duplicate reads
            self.assertEqual(int(run_rows[0]["N_LANES"]), 200)

            with (root / "spread_out.read_stats.tsv").open() as handle:
                summary = list(csv.DictReader(handle, delimiter="\t"))[0]
            mean_gc = float(summary["MEAN_READ_GC_PERC"])
            self.assertGreaterEqual(mean_gc, 0.35)
            self.assertLessEqual(mean_gc, 0.65)

    def test_bam_multiwindow_high_density_and_sparse_tail_backfill(self):
        with tempfile.TemporaryDirectory(prefix="qrs-bam-edge-") as tmp:
            root = Path(tmp)
            sam = root / "edge.sam"
            bam = root / "edge.bam"
            prefix = root / "edge_out"

            # chr1: LN=100 with 500 reads (mapped > chrom_len, multiple reads per coordinate)
            # chr2: LN=1000000 with 500 reads clustered in first 500 bp (99.95% empty tail -> backfill)
            header = [
                "@HD\tVN:1.6\tSO:coordinate",
                "@SQ\tSN:chr1\tLN:100",
                "@SQ\tSN:chr2\tLN:1000000",
                "@RG\tID:rg1\tSM:edge_sample",
            ]
            records = []
            uid = 0
            for pos in range(1, 91):
                for _ in range(5):
                    records.append(
                        f"A01:1:FC1:U{uid}:1:1:1\t0\tchr1\t{pos}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1\tNM:i:0"
                    )
                    uid += 1
            for pos in range(1, 501):
                records.append(
                    f"A01:1:FC1:U{uid}:1:1:1\t0\tchr2\t{pos}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1\tNM:i:0"
                )
                uid += 1

            sam.write_text("\n".join(header + records) + "\n")
            subprocess.run([samtools_bin(), "view", "-b", "-o", str(bam), str(sam)],
                           check=True, capture_output=True)
            subprocess.run([samtools_bin(), "index", str(bam)],
                           check=True, capture_output=True)

            run = subprocess.run(
                [qrs_bin(), "--bam", str(bam), "--subset", "300",
                 "--write_read_gc", "--write_run_stats", "--out", str(prefix)],
                text=True, capture_output=True,
            )
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)

            with (root / "edge_out.GC_stats/edge_sample.GC_stats.tsv").open() as handle:
                gc_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(gc_rows), 300)

            with (root / "edge_out.per_sample_run.tsv").open() as handle:
                run_rows = list(csv.DictReader(handle, delimiter="\t"))
            # Proves all 300 reads across Pass 1 + Pass 2 backfill are distinct
            self.assertEqual(int(run_rows[0]["N_LANES"]), 300)


if __name__ == "__main__":
    unittest.main()

