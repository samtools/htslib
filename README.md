HTSlib is an experimental implementation of a unified C library for accessing
common file formats, such as [SAM][1] and [VCF][2], used for high-throughput
sequencing data.
HTSlib only depends on [zlib][3].
It is known to be compatible with gcc, g++ and clang.

HTSlib implements a generalized BAM index, with
file extension `.csi` (coordinate-sorted index). The HTSlib file reader first
looks for the new index and then for the old if the new index is absent.

HTSlib is unfinished. It has not been tested on large-scale real data. Some
useful APIs are missing.

[1]: http://samtools.sourceforge.net
[2]: http://vcftools.sourceforge.net/specs.html
[3]: http://zlib.net/
