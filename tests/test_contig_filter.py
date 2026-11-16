"""Run with QRS_BIN=/path/to/qrs SAMTOOLS=/path/to/samtools python3 -m unittest discover -s tests -v."""

import csv
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


class ContigFilterTest(unittest.TestCase):
    def check_contigs(self, contigs, expected_reads, subset=-1, write_run_stats=False):
        with tempfile.TemporaryDirectory(prefix="qrs-contigs-") as tmp:
            root = Path(tmp)
            sam = root / "input.sam"
            bam = root / "input.bam"
            prefix = root / "stats"
            header = ["@HD\tVN:1.6\tSO:coordinate"]
            header.extend(f"@SQ\tSN:{name}\tLN:1000" for name, _ in contigs)
            header.append("@RG\tID:rg1\tSM:sample1")
            records = []
            for idx, (name, canonical) in enumerate(contigs):
                sequence = "ACGT" if canonical else "GGCC"
                for pos in (2, 10):
                    records.append(
                        f"INST:1:FLOW:1:{idx}:{pos}:1\t0\t{name}\t{pos}\t60"
                        f"\t4M\t*\t0\t0\t{sequence}\tIIII\tRG:Z:rg1"
                    )
            sam.write_text("\n".join(header + records) + "\n")
            samtools = os.environ.get("SAMTOOLS", "samtools")
            subprocess.run([samtools, "view", "-b", "-o", str(bam), str(sam)],
                           check=True, capture_output=True)
            subprocess.run([samtools, "index", str(bam)], check=True, capture_output=True)
            cmd = [
                os.environ.get("QRS_BIN", "qrs"),
                f"--subset={subset}",
                "--write_read_gc",
                "--out", str(prefix),
                "--bam", str(bam),
            ]
            if write_run_stats:
                cmd.append("--write_run_stats")
            run = subprocess.run(cmd, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertIn(f"Mapped reads: {expected_reads}", run.stdout)
            with (root / "stats.GC_stats/sample1.GC_stats.tsv").open() as handle:
                gc_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(gc_rows), expected_reads)
            self.assertTrue(all(float(row["GC_PERC"]) == 0.5 for row in gc_rows))
            with (root / "stats.read_stats.tsv").open() as handle:
                summary = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(summary), int(expected_reads > 0))
            if summary:
                self.assertEqual(float(summary[0]["MEAN_READ_GC_PERC"]), 0.5)

            self.assertEqual((root / "stats.rundetails.tsv").exists(), write_run_stats)
            self.assertEqual((root / "stats.per_sample_run.tsv").exists(), write_run_stats)

    def test_mixed_contigs(self):
        # Consecutive exclusions and retained contigs with nonconsecutive BAM IDs.
        contigs = [("chrM", False), ("chrUn_a", False), ("chr1", True),
                   ("chr1_alt", False), ("chrX", True), ("2", True),
                   ("GL0001", False)]
        for subset in (-1, 6):
            with self.subTest(subset=subset):
                self.check_contigs(contigs, expected_reads=6, subset=subset)

    def test_all_contigs_excluded(self):
        self.check_contigs([("chrM", False), ("chrUn_a", False)], expected_reads=0)

    def test_all_contigs_retained(self):
        self.check_contigs([("chr1", True), ("chrX", True), ("2", True)],
                           expected_reads=6, write_run_stats=True)


if __name__ == "__main__":
    unittest.main()
