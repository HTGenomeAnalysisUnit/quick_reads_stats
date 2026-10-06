# Package

version       = "0.1.1"
author        = "edoardo.giacopuzzi"
description   = "Quickly collect essential read-level and alignment stats from Illumina NGS BAM files"
license       = "MIT"
srcDir        = "src"
bin           = @["qrs"]
skipDirs      = @["test"]


# Dependencies

requires "nim >= 1.4.8", "hts >= 0.3.21", "argparse >= 3.0.0"
