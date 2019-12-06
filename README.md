[![Build Status](https://travis-ci.org/samtools/htslib.svg?branch=develop)](https://travis-ci.org/samtools/htslib)
[![Build status](https://ci.appveyor.com/api/projects/status/v46hkwyfjp3l8nd3/branch/develop?svg=true)](https://ci.appveyor.com/project/samtools/htslib/branch/develop)
[![Github All Releases](https://img.shields.io/github/downloads/samtools/htslib/total.svg)](https://github.com/samtools/htslib)

HTSlib is an implementation of a unified C library for accessing common file
formats, such as [SAM, CRAM and VCF][1], used for high-throughput sequencing
data, and is the core library used by [samtools][2] and [bcftools][3].
HTSlib only depends on [zlib][4].
It is known to be compatible with gcc, g++ and clang.

HTSlib implements a generalized BAM index, with file extension `.csi`
(coordinate-sorted index). The HTSlib file reader first looks for the new index
and then for the old if the new index is absent.

This project also includes the popular tabix indexer, which indexes both `.tbi`
and `.csi` formats, and the bgzip compression utility.

[1]: http://samtools.github.io/hts-specs/
[2]: http://github.com/samtools/samtools
[3]: http://samtools.github.io/bcftools/
[4]: http://zlib.net/

### Building HTSlib

See [INSTALL](INSTALL) for complete details.
[Release tarballs][download] contain generated files that have not been
committed to this repository, so building the code from a Git repository
requires extra steps:

```sh
autoheader     # If using configure, generate the header template...
autoconf       # ...and configure script (or use autoreconf to do both)
./configure    # Optional but recommended, for choosing extra functionality
make
make install
```

[download]: http://www.htslib.org/download/
