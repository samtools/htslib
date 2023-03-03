[![Github All Releases](https://img.shields.io/github/downloads/samtools/htslib/total.svg)](https://github.com/samtools/htslib/samples)

HTSlib is an implementation of a unified C library for accessing common file
formats, such as [SAM, CRAM and VCF][1], used for high-throughput sequencing
data, and is the core library used by [samtools][2] and [bcftools][3].

A set of sample programs are available which showcases the usage of APIs in HTSlib.
They are based on version 1.17 of HTSLib and are mainly for demonstration of API usage.
Further optimization and error handling might be required for actual usage.


[1]: http://samtools.github.io/hts-specs/
[2]: http://github.com/samtools/samtools
[3]: http://samtools.github.io/bcftools/

### Building and using sample programs

GCC and compatible compilers can be used to build these samples.

A makefile is available along with source files which links statically to
htslib. To use dynamic linking, update the makefile's 'LDFLAGS' and 'rpath'
path. The 'rpath' path to be set as the path to lib directory of htslib
installation.

```sh

# linking statically on a linux machine
gcc -g -o <binary name> <source name> -I <htslib include path> \
    <htslib lib path>/libhts.a -lcrypto -lm -lpthread -lcurl -llzma -lz -lbz2

# dynamically linking with libhts
gcc -g -o <binary name> <source name> -I <htslib include path> \
    -L <htslib lib path> -lhts -Wl,-rpath,<htslib lib path>

```

In many cases, the alignment data are expected as sorted, compressed and
indexed.

### The samples...

[Flags][Flags]

  This application showcases the basic read of alignment files and flag
  access. It reads and shows the count of read1 and read2 alignments.

[Split][Split]

  This application showcases the basic read and write of alignment data. It
  saves the read1 and read2 as separate files in given directory, one as sam
  and other as bam.

[Split2][Split2]

  This application showcases the output file format selection. It saves the
  read1 and read2 as separate files in given directory, both as compressed
  sam though the extensions are different.

[Cram][Cram]

  This application showcases the different way in which cram reference data
  is used for cram output creation.

[Read_fast][Read_fast]

  This application showcases the fasta/fastq data read.

[Read_header][Read_header]

  This application showcases the read and access of header data. It can show
  all header line of given type, data of a given tag on a specific header
  line or for all lines of given type.

[Read_ref][Read_ref]

  This application showcases the read and access of header data. It shows
  all reference names which has length equal or greather to given input.

[Read_bam][Read_bam]

  This application showcases read of different alignment data fields. It
  shows contents of each alignment.

[Read_aux][Read_aux]

  This application showcases read of specific auxiliary tag data in
  alignment. It shows the data retrieved using 2 APIs, one as a string with
  tag data and other as raw data alternatively.

[Dump_aux][Dump_aux]

  This application showcases read of all auxiliary tag data one by one in an
  alignment. It shows the data retrieved.

[Add_header][Add_header]

  This application showcases the write of header lines to a file. It adds
  header line of types, SQ, RG, PG and CO and writes to standard output.

[Remove_header][Remove_header]

  This application showcases removal of header line from a file. It removes
  either all header lines of given type or one specific line of given type
  with given unique identifier. Modified header is written on standard
  output.

[Update_header][Update_header]

  This application shows the update of header line fields, where update is
  allowed. It takes the header line type, unique identifier for the line,
  tag to be modified and the new value. Updated data is written on standard
  output.

[Mod_bam][Mod_bam]

  This application showcases the update of alignment data. It takes
  alignment name, position of field to be modified and new value of
  it. Modified data is written on standard output.

[Mod_aux][Mod_aux]

  This application showcases the update of auxiliary data in alignment. It
  takes alignment name, tag to be modified, its type and new value. Modified
  data is written on standard output.

[Mod_aux_ba][Mod_aux_ba]

  This application showcases the update of auxiliary array data in
  alignment. It adds count of ATCGN base as an array in auxiliary data,
  BA:I. Modified data is written on standard output.

[Write_fast][Write_fast]

  This application showcases the fasta/fastq data write. It appends a dummy
  data to given file.

[Index_write][Index_write]

  This application showcases the creation of index along with output
  creation. Based on file type and shift, it creates bai, csi or crai files.

[Read_reg][Read_reg]:

  This application showcases the usage of region specification in alignment
  read.

[Read_multireg][Read_multireg]:

  This application showcases the usage of mulitple region specification in
  alignment read.

[Pileup][Pileup]:

  This application showcases the pileup api, where all alignments covering a
  reference position are accessed together. It displays the bases covering
  each position on standard output.

[Mpileup][Mpileup]:

  This application showcases the mpileup api, which supports multiple input
  files for pileup and gives a side by side view of them in pileup
  format. It displays the bases covering each position on standard output.

[Modstate][Modstate]:

  This application showcases the access of base modifications in
  alignment. It shows the modifications present in an alignment and accesses
  them using available APIs. There are 2 APIs and which one to be used can
  be selected through input.

[Pileup_mod][Pileup_mod]:

  This application showcases the base modification access in pileup mode. It
  shows the pileup display with base modifications.

[Flags_field][Flags_field]

  This application showcases the read of selected fields alone, reducing the
  overhead / increasing the performance. It reads the flag field alone and
  shows the count of read1 and read2. This has impact only on CRAM files.

[Split_thread1][Split_thread1]

  This application showcases the use of threads in file handling. It saves
  the read1 and read2 as separate files in given directory, one as sam and
  other as bam. 2 threads are used for read and 1 each dedicated for each
  output file.

[Split_thread2][Split_thread2]

  This application showcases the use of thread pool in file handling. It
  saves the read1 and read2 as separate files in given directory, one as sam
  and other as bam. A pool of 4 threads is created and shared for both read
  and write.

### More Information

More detailed documentation is available in the [DEMO.md][DEMO] with worked
examples per demonstration tool.


[Flags]: flags_demo.c
[Split]: split.c
[Split2]: split2.c
[Cram]: cram.c
[Read_fast]: read_fast.c
[Read_header]: read_header.c
[Read_ref]: read_refname.c
[Read_bam]: read_bam.c
[Read_aux]: read_aux.c
[Dump_aux]: dump_aux.c
[Add_header]: add_header.c
[Remove_header]: rem_header.c
[Update_header]: update_header.c
[Mod_bam]: mod_bam.c
[Mod_aux]: mod_aux.c
[Mod_aux_ba]: mod_aux_ba.c
[Write_fast]: write_fast.c
[Index_write]: index_write.c
[Read_reg]: index_reg_read.c
[Read_multireg]: index_multireg_read.c
[Pileup]: pileup.c
[Mpileup]: mpileup.c
[Modstate]: modstate.c
[Pileup_mod]: pileup_mod.c
[Flags_field]: flags_htsopt_field.c
[Split_thread1]: split_thread1.c
[Split_thread2]: split_thread2.c
[DEMO]: DEMO.md
