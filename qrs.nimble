# Package

version       = "0.1.0"
author        = "edoardo.giacopuzzi"
description   = "Illumina run stats from read names"
license       = "MIT"
srcDir        = "src"
bin           = @["qrs"]
skipDirs      = @["test"]


# Dependencies

requires "nim >= 1.4.8", "hts >= 0.3.21", "argparse >= 3.0.0"
