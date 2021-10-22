[![Build Status](https://api.cirrus-ci.com/github/samtools/htslib.svg?branch=develop)](https://api.cirrus-ci.com/github/samtools/htslib)
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

This project also includes the popular tabix indexer, which creates both `.tbi`
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
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
make install
```

[download]: http://www.htslib.org/download/

### Citing

Please cite this paper when using HTSlib for your publications.

> HTSlib: C library for reading/writing high-throughput sequencing data </br>
> James K Bonfield, John Marshall, Petr Danecek, Heng Li, Valeriu Ohan, Andrew Whitwham, Thomas Keane, Robert M Davies </br>
> _GigaScience_, Volume 10, Issue 2, February 2021, giab007, https://doi.org/10.1093/gigascience/giab007

```
@article{10.1093/gigascience/giab007,
    author = {Bonfield, James K and Marshall, John and Danecek, Petr and Li, Heng and Ohan, Valeriu and Whitwham, Andrew and Keane, Thomas and Davies, Robert M},
    title = "{HTSlib: C library for reading/writing high-throughput sequencing data}",
    journal = {GigaScience},
    volume = {10},
    number = {2},
    year = {2021},
    month = {02},
    abstract = "{Since the original publication of the VCF and SAM formats, an explosion of software tools have been created to process these data files. To facilitate this a library was produced out of the original SAMtools implementation, with a focus on performance and robustness. The file formats themselves have become international standards under the jurisdiction of the Global Alliance for Genomics and Health.We present a software library for providing programmatic access to sequencing alignment and variant formats. It was born out of the widely used SAMtools and BCFtools applications. Considerable improvements have been made to the original code plus many new features including newer access protocols, the addition of the CRAM file format, better indexing and iterators, and better use of threading.Since the original Samtools release, performance has been considerably improved, with a BAM read-write loop running 5 times faster and BAM to SAM conversion 13 times faster (both using 16 threads, compared to Samtools 0.1.19). Widespread adoption has seen HTSlib downloaded \\&gt;1 million times from GitHub and conda. The C library has been used directly by an estimated 900 GitHub projects and has been incorporated into Perl, Python, Rust, and R, significantly expanding the number of uses via other languages. HTSlib is open source and is freely available from htslib.org under MIT/BSD license.}",
    issn = {2047-217X},
    doi = {10.1093/gigascience/giab007},
    url = {https://doi.org/10.1093/gigascience/giab007},
    note = {giab007},
    eprint = {https://academic.oup.com/gigascience/article-pdf/10/2/giab007/36332285/giab007.pdf},
}
```
