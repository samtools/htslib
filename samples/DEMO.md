# HTS API

## HTSLib APIs and samtools

HTSLib is a C library implementation used to access and process the genome
sequence data. HTSLib implements multiple API interfaces, HTS API, VCF API and
SAM API.  HTS API provides a framework for use by other APIs and applications,
implements bgzf compression, htscodecs and provides CRAM format support.  VCF
APIs work with variant data in VCF and BCF format.

SAM API works with sequence data of different formats, SAM / BAM / CRAM /
FASTA / FASTQ, and provides methods to do operations on the data. It uses
methods from HTS API.

'samtools' is the utility used to read and modify sequence data. It uses SAM
APIs from HTSLib to work on the sequence data.


## About this document

There are a number of demonstration utilities and their source code in
'samples' directory of HTSLib and this document gives the description of them
and the usage of API of HTSLib.  The samples are for demonstration
purposes only and proper error handling is required for actual usage. This
document is based on HTSLib version 1.17.

Updates to this document may be made along with later releases when required.


## The sample apps

Flags - This application showcases the basic read of alignment files and flag
access. It reads and shows the count of read1 and read2 alignments.

Split - This application showcases the basic read and write of alignment data.
It saves the read1 and read2 as separate files in given directory, one as sam
and other as bam.

Split2 - This application showcases the output file format selection. It saves
the read1 and read2 as separate files in given directory, both as compressed
sam though the extensions are different.

Cram - This application showcases the different way in which cram reference
data is used for cram output creation.

Read_fast - This application showcases the fasta/fastq data read.

Read_header - This application showcases the read and access of header data.
It can show all header line of given type, data of a given tag on a specific
header line or for all lines of given type.

Read_ref - This application showcases the read and access of header data.
It shows all reference names which has length equal or greater to given input.

Read_bam - This application showcases read of different alignment data fields.
It shows contents of each alignment.

Read_aux - This application showcases read of specific auxiliary tag data in
alignment. It shows the data retrieved using 2 APIs, one as a string with tag
data and other as raw data alternatively.

Dump_aux - This application showcases read of all auxiliary tag data one by one
in an alignment. It shows the data retrieved.

Add_header - This application showcases the write of header lines to a file.
It adds header line of types, SQ, RG, PG and CO and writes to standard output.

Remove_header - This application showcases removal of header line from a file.
It removes either all header lines of given type or one specific line of given
type with given unique identifier. Modified header is written on standard
output.

Update_header - This application shows the update of header line fields, where
update is allowed. It takes the header line type, unique identifier for the
line, tag to be modified and the new value. Updated data is written on standard
output.

Mod_bam - This application showcases the update of alignment data. It takes
alignment name, position of field to be modified and new value of it.
Modified data is written on standard output.

Mod_aux - This application showcases the update of auxiliary data in alignment.
It takes alignment name, tag to be modified, its type and new value. Modified
data is written on standard output.

Mod_aux_ba - This application showcases the update of auxiliary array data in
alignment. It adds count of ATCGN base as an array in auxiliary data, BA:I.
Modified data is written on standard output.

Write_fast - This application showcases the fasta/fastq data write. It appends
data to given file.

Index_write - This application showcases the creation of index along with
output creation. Based on file type and shift, it creates bai, csi or crai
files.

Index_fast - This application showcases the index creation on fasta/fastq
reference files.

Read_reg - This application showcases the usage of region specification in
alignment read.

Read_multireg - This application showcases the usage of multiple region
specification in alignment read.

Read_fast_index - This application showcases the fasta/fastq data read using
index.

Pileup - This application showcases the pileup api, where all alignments
covering a reference position are accessed together. It displays the bases
covering each position on standard output.

Mpileup - This application showcases the mpileup api, which supports multiple
input files for pileup and gives a side by side view of them in pileup format.
It displays the bases covering each position on standard output.

Modstate - This application showcases the access of base modifications in
alignment. It shows the modifications present in an alignment and accesses them
using available APIs. There are 2 APIs and which one to be used can be selected
through input.

Pileup_mod - This application showcases the base modification access in pileup
mode. It shows the pileup display with base modifications.

Flags_field - This application showcases the read of selected fields alone,
reducing the overhead / increasing the performance. It reads the flag field
alone and shows the count of read1 and read2. This has impact only on CRAM
files.

Split_thread1 - This application showcases the use of threads in file handling.
It saves the read1 and read2 as separate files in given directory, one as sam
and other as bam. 2 threads are used for read and 1 each dedicated for each
output file.

Split_thread2 - This application showcases the use of thread pool in file
handling. It saves the read1 and read2 as separate files in given directory,
one as sam and other as bam. A pool of 4 threads is created and shared for both
read and write.

Qtask_ordered - This application showcases the use of queues and threads for
custom processing. Alignments in input file are updated with their GC ratio
on a custom aux tag. The processing may occur in any order but the result is
retrieved in same order as it was queued and saved to disk.

Qtask_unordered - This application showcases the use of queues and threads
for custom processing. The count of bases and GC ratio are calculated and
displayed.  The order of counting is irrelevant and hence ordered retrieval is
not used.

## Building the sample apps

The samples expect the HTSLib is installed, libraries and header file path are
part of the PATH environment variable. If not, these paths need to be explicitly
passed during the build time.

Gcc and compatible compilers can be used to build the samples.

These applications can be linked statically or dynamically to HTSLib.
For static linking, along with htslib other libraries and/or headers required
to build are, math, pthread, curl, lzma, z and bz2 libraries.

A makefile is available along with source files which links statically to
htslib. To use dynamic linking, update the makefile's 'LDFLAGS' and 'rpath'
path. The 'rpath' path to be set as the path to lib directory of htslib
installation.


## Usage of HTS APIs
### Sequence data file access for read

The sequence data file for read may be opened using the sam_open method. It
opens the file and returns samFile (htsFile) pointer on success or NULL on
failure.  The input can be path to a file in disk, network, cloud or '-'
designating the standard input.

SAM, BAM and CRAM file formats are supported and the input file format is
detected from the file content.

Once done with the file, it needs to be closed with sam_close.

Many times, header details would be required and can be read using
sam_hdr_read api.  It returns sam_hdr_t pointer or NULL. The returned header
needs to be destroyed using sam_hdr_destroy when no longer required.

The sequence data may be compressed or uncompressed on disk and on memory it
is read and kept as uncompressed BAM format. It can be read from a file using
sam_read1 api. samFile pointer, header and bam storage are to be passed as
argument and it returns 0 on success, -1 on end of file and < -1 in case of
errors.

The bam storage has to be initialized using bam_init1 api before the call and
can be reused for successive reads. Once done, it needs to be destroyed using
bam_destroy1.  The member field named core - bam1_core_t - in bam storage,
bam1_t, has the sequence data in an easily accessible way. Using the fields
and macros, data can easily be read from it.

    #include <htslib/sam.h>

    int main(int argc, char *argv[])
    {
        ...
        //initialize
        if (!(bamdata = bam_init1()))
           ... // error
        //open input files - r reading
        if (!(infile = sam_open(inname, "r")))
           ... // error
        //read header
        if (!(in_samhdr = sam_hdr_read(infile)))
           ... // error

        //read data, check flags and update count
        while ((c = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
            if (bamdata->core.flag & BAM_FREAD1)
                cntread1++;
            ...

        //clean up
        if (in_samhdr)
            sam_hdr_destroy(in_samhdr);

        if (infile)
            sam_close(infile);

        if (bamdata)
            bam_destroy1(bamdata);

        return ret;
    }
Refer: flags_demo.c

This shows the count of read1 and read2 alignments.

    ./flags /tmp/sample.sam.gz

To read CRAM files, reference data is required and if it is not available, based
on configuration, library may try to download it from external repositories.


### Sequence data file access for write

File access for write is similar to read with a few additional optional steps.

The output file can be opened using sam_open api as in read, with "w" instead
of "r" as mode. This opens the file for writing and uses mode to select the
output file type. "w" alone denotes SAM, "wb" denotes BAM and "wc" denotes CRAM.

Another way is to use sam_open_mode method, which sets the output file type and
compression based on the file name and explicit textual format specification.
This method expects a buffer to append type and compression flags. Usually a
buffer with standard file open flag is used, the buffer past the flag is passed
to the method to ensure existing flags and updates from this method are present
in the same buffer without being overwritten. This method will add more flags
indicating file type and compression based on name. If explicit format detail
given, then extension is ignored and the explicit specification is used. This
updated buffer can be used with sam_open to select the file format.

sam_open_format method may also be used to open the file for output as more
information on the output file can be specified using this. Can use
mode buffer from sam_open_mode api or explicit format structure for this.

The header data can be written using the sam_hdr_write api. When the header
data is copied to another variable and has different lifetime, it is good to
increase the reference count of the header using sam_hdr_incr_ref and
sam_hdr_destroy called as many times as required.

The alignment data can be written using the sam_write1 api. It takes a samFile
pointer, header pointer and the alignment data. The header data is required to
set the reference name in the alignment. It returns -ve value on error.

    int main(int argc, char *argv[])
    {
        ...
        if (!(infile = sam_open(inname, "r")))
           ... // error
        outfile1 = sam_open(file1, "w");            //as SAM
        outfile2 = sam_open(file2, "wb");           //as BAM
        ...
        if (!(in_samhdr = sam_hdr_read(infile)))
           ... // error

        //write header
        if ((sam_hdr_write(outfile1, in_samhdr) == -1) ||
         (sam_hdr_write(outfile2, in_samhdr) == -1))
           ... // error

        while ((c = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
            if (bamdata->core.flag & BAM_FREAD1) {
                if (sam_write1(outfile1, in_samhdr, bamdata) < 0) {
                    ... // error
    }
Refer: split.c

This creates 1.sam and 2.bam in /tmp/ containing read1 and read2 respectively.

    ./split /tmp/sample.sam.gz /tmp/

Below code excerpt shows sam_open_mode api usage.

    int main(int argc, char *argv[])
    {
        ...
        //set file open mode based on file name for 1st and as explicit for 2nd
        if ((sam_open_mode(mode1+1, file1, NULL) == -1) ||
         (sam_open_mode(mode2+1, file2, "sam.gz") == -1))
           ... // error
        if (!(infile = sam_open(inname, "r")))
           ... // error

        //open output files
        outfile1 = sam_open(file1, mode1);                          //as compressed SAM through sam_open
        outfile2 = sam_open_format(file2, mode2, NULL);             //as compressed SAM through sam_open_format
        ...
    }
Refer: split2.c

This creates 1.sam.gz and 2.sam in /tmp/ both having compressed data.

    ./split2 /tmp/sample.sam.gz /tmp/

An htsFormat structure filled appropriately can also be used to specify output
file format while using sam_open_format api.


### CRAM writing

CRAM files uses reference data and compresses alignment data. A CRAM file may
be created with external reference data file - most appropriate, with embedded
reference in it or with no reference data at all. It can also be created using
an autogenerated reference, based on consensus with-in the alignment data.
The reference detail can be set to an htsFormat structure using hts_parse_format
api and used with sam_open_format api to create appropriate CRAM file.

    ...
    snprintf(reffmt1, size1, "cram,reference=%s", reffile);
    snprintf(reffmt2, size2, "cram,embed_ref=1,reference=%s", reffile);
    ...
    if (hts_parse_format(&fmt1, reffmt1) == -1 ||                   //using external reference - uses the M5/UR tags to get
     reference data during read
            hts_parse_format(&fmt2, reffmt2) == -1 ||               //embed the reference internally
            hts_parse_format(&fmt3, "cram,embed_ref=2") == -1 ||    //embed autogenerated reference
            hts_parse_format(&fmt4, "cram,no_ref=1") == -1) {       //no reference data encoding at all
       ... // error
    outfile1 = sam_open_format(file1, "wc", &fmt1); outfile2 = sam_open_format(file2, "wc", &fmt2);
    ...
Refer: cram.c


### FASTA/FASTQ data access

FASTA/FASTQ files have the raw sequence data and the data can be read one by
one using sam_read1 or a selected range using a region. The data can be written
similar to alignment data using sam_write1 api. To write the file, format
can be set by updating mode buffer using sam_open_mode with file name
or explicit format text. This mode buffer can be used with sam_open or can be
used with sam_open_format with explicit format information in htsFormat
structure.

It is the FASTA format which is mainly in use to store the reference data.

    ...
    if (!(bamdata = bam_init1()))
      ... // error
    if (!(infile = sam_open(inname, "r")))
       ... // error
    if (infile->format.format != fasta_format && infile->format.format != fastq_format)
       ... // error
    if (!(in_samhdr = sam_hdr_read(infile)))
       ... // error

    while ((c = sam_read1(infile, in_samhdr, bamdata)) >= 0)
       ... // error
        printf("\nsequence: ");
        for (c = 0; c < bamdata->core.l_qseq; ++c) {
            printf("%c", seq_nt16_str[bam_seqi(bam_get_seq(bamdata), c)]);
        }
        if (infile->format.format == fastq_format) {
            printf("\nquality: ");
            for (c = 0; c < bamdata->core.l_qseq; ++c) {
                printf("%c", bam_get_qual(bamdata)[c] + 33);
    ...
Refer: read_fast.c

    ...
    char mode[4] = "a";
    ...
    if (sam_open_mode(mode + 1, outname, NULL) < 0)
       ... // error
    if (!(outfile = sam_open(outname, mode)))
       ... // error
    if (bam_set1(bamdata, strlen(name), name, BAM_FUNMAP, -1, -1, 0, 0, NULL, -1, -1, 0, strlen(data), data, qual, 0) < 0)
       ... // error
    if (sam_write1(outfile, out_samhdr, bamdata) < 0) {
        printf("Failed to write data\n");
        ...
Refer: write_fast.c


### Header data read

The header gives the version, reference details, read group, change history
and comments. These data are stored inside the sam_hdr_t. Each of these
entries, except comments, have their unique identifier and it is required to
access different fields of them.  The api sam_hdr_count_lines gives the count
of the specified type of header line.  The value of a unique identifier to a
specific type of header line can be retrieved with sam_hdr_line_name api.  The
api sam_hdr_find_tag_id and sam_hdr_find_tag_pos can get the field data from a
header line using unique identifier values or using position. The full header
line can be retrieved using sam_hdr_find_line_pos or sam_hdr_line_id with
position and unique identifier values respectively.

    ...
    if (!(in_samhdr = sam_hdr_read(infile)))
        ... // error
    ...
      if (tag)
          ret = sam_hdr_find_tag_id(in_samhdr, header, id, idval, tag, &data);
      else
          ret = sam_hdr_find_line_id(in_samhdr, header, id, idval, &data);
    ...
        linecnt = sam_hdr_count_lines(in_samhdr, header);
        ...
            if (tag)
                ret = sam_hdr_find_tag_pos(in_samhdr, header, c, tag, &data);
            else
                ret = sam_hdr_find_line_pos(in_samhdr, header, c, &data);
        ...
Refer: read_header.c

This will show the VN tag's value from HD header.

    ./read_header /tmp/sample.sam.gz HD VN

Shows the 2nd SQ line's LN field value.

    ./read_header /tmp/sample.sam.gz SQ SN T2 LN

Below code excerpt shows the reference names which has length above given value.

    ...
    linecnt = sam_hdr_count_lines(in_samhdr, "SQ"); //get reference count
    ...
    //iterate and check each reference's length
    for (pos = 1, c = 0; c < linecnt; ++c) {
        if ((ret = sam_hdr_find_tag_pos(in_samhdr, "SQ", c, "LN", &data) == -2))
            ... // error

        size = atoll(data.s);
        if (size < minsize) {
            //not required
            continue;
        }

        //sam_hdr_find_tag_pos(in_samhdr, "SQ", c, "SN", &data) can also do the same!
        if (!(id = sam_hdr_line_name(in_samhdr, "SQ", c)))
            ... // error

        printf("%d,%s,%s\n", pos, id, data.s);
    ...
Refer: read_refname.c


### Alignment data read

The alignment / sequence data contains many fields. Mainly the read/query
name, flags indicating the properties of the read, reference sequence name,
position in reference to which it matches, quality of the read, CIGAR string
indicating the match status, position of mate / reverse strand, name of
reference sequence to which mate matches, the insert length, base sequence,
quality value of each base and auxiliary fields.

Header data would be required to retrieve the reference names as alignment
contains the position of the reference in the header.

A few of the data are directly visible in bam1_t and the rest are hidden
inside data member of bam1_t and can easily be retrieved using macros.
bam_get_qname gives the name of the read, sam_hdr_tid2name gives the reference
name. bam_get_cigar retrieves the cigar operation array, which can be decoded
using bam_cigar_oplen to get count of bases to which that operation applicable
and bam_cigar_opchr to get the cigar operation.  bam_seqi retrieves the base
data at a given position in alignment and it can be converted to character by
indexing the seq_nt16_str array.

    ...
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0)
    {
        //QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAG:TYPE:VALUE]
        printf("NAME: %s\n", bam_get_qname(bamdata));                           //get the query name using the macro
        flags = bam_flag2str(bamdata->core.flag);                               //flags as string
        ...
        tidname = sam_hdr_tid2name(in_samhdr, bamdata->core.tid);
        ...
        printf("MQUAL: %d\n", bamdata->core.qual);                              //map quality value
        cigar = bam_get_cigar(bamdata);                                         //retrieves the cigar data
        for (i = 0; i < bamdata->core.n_cigar; ++i) {                           //no. of cigar data entries
            printf("%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
            //the macros gives the count of operation and the symbol of operation for given cigar entry
        }
        printf("\nTLEN/ISIZE: %"PRIhts_pos"\n", bamdata->core.isize);
        data = bam_get_seq(bamdata);
                 //get the sequence data
        if (bamdata->core.l_qseq != bam_cigar2qlen(bamdata->core.n_cigar, cigar)) { //checks the length with CIGAR and query
        ...
        for (i = 0; i < bamdata->core.l_qseq ; ++i) {       //sequence length
            printf("%c", seq_nt16_str[bam_seqi(data, i)]);  //retrieves the base from (internal compressed) sequence data
            ...
            printf("%c", bam_get_qual(bamdata)[i]+33);      //retrieves the quality value
        ...
Refer: read_bam.c

Shows the data from alignments.

    ./read_bam /tmp/sample.sam.gz


### Aux data read

Auxiliary data gives extra information about the alignment. There can be a
number of such data and can be accessed by specifying required tag or by
iterating one by one through them once the alignment is read as bam1_t. The
auxiliary data are stored along with the variable length data in the data
field of bam1_t. There are macros defined to retrieve information about
auxiliary data from the data field of bam1_t.

Data for a specific tag can be retrieved as a string or can be retrieved as raw
data. bam_aux_get_str retrieves as a string, with tag name, tag type and data.
bam_aux_get can get raw data and with bam_aux_type and bam_aux2A, bam_aux2f etc.
the raw data can be extracted.

To iterate through all data, the start of aux data is retrieved using macro
bam_aux_first and successive ones using bam_aux_next. Macro bam_aux_tag gives
the tag of the aux field and bam_aux_type gives the information about type of
the aux field.

Bam_aux2i, bam_aux2f, bam_aux2Z macros retrieve the aux data's value as
integer, float and string respectively. The integer value may be of different
precision / size and the bam_aux_type character indicates how to use the
value. The string/hex data are NULL terminated.

For array data, bam_aux_type will return 'B' and bam_auxB_len gives the length
of the array. bam_aux_type with the next byte will give the type of data in
the array. bam_auxB2i, bam_auxB2f will give integer and float data from a
given position of the array.

    ...
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        //option 1 - get data as string with tag and type
        if ((c = bam_aux_get_str(bamdata, tag, &sdata)) == 1) {
            printf("%s\n",sdata.s);
        ...
        //option 2 - get raw data
        if ((data = bam_aux_get(bamdata, tag)) != NULL) {
            printauxdata(stdout, bam_aux_type(data), -1, data);
    ...
Refer: read_aux.c

Shows the MD aux tag from alignments.

    ./read_aux ../../samtools/test/mpileup/mpileup.1.bam MD

    ...
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        data = bam_aux_first(bamdata);                                              //get the first aux data
        while (data) {
            printf("%.2s:%c:", bam_aux_tag(data), NULL != strchr("cCsSiI", bam_aux_type(data)) ? 'i' : bam_aux_type(data));
              //macros gets the tag and type of aux data
            //dump the data
            printauxdata(stdout, bam_aux_type(data), -1, data);
            ...
            data = bam_aux_next(bamdata, data);                                     //get the next aux data
    ...
Refer: dump_aux.c

Shows all the tags from all alignments.

    ./dump_aux ../../samtools/test/mpileup/mpileup.1.bam


### Add/Remove/Update header

There are specific types of data that can be part of header data. They have
a tag from HD, SQ, RG, PG and CO. Fully formatted header lines, separated by new
line, can be added with sam_hdr_add_lines api. A single header line can be added
using sam_hdr_add_line api where the header type, tag and value pair are passed
as arguments, terminated with a NULL argument. The PG header lines are special
that they have a kind of linkage to previous PG lines. This linkage can be auto
generated by using sam_hdr_add_pg api which sets the 'PP' field used in linkage.
sam_hdr_write api does the write of the header data to file.

    ...
    //add SQ line with SN as TR1 and TR2
    if (sam_hdr_add_lines(in_samhdr, &sq[0], 0))                                        //length as 0 for NULL terminated data
        ... // error

    //add RG line with ID as RG1
    if (sam_hdr_add_line(in_samhdr, "RG", "ID", "RG1", "LB", "Test", "SM", "S1", NULL))
        ... // error

    //add PG/CO lines
    if (sam_hdr_add_pg(in_samhdr, "add_header", "VN", "Test", "CL", data.s, NULL))      //NULL is to indicate end of args
        ... // error
    if (sam_hdr_add_line(in_samhdr, "CO", "Test data", NULL))                           //NULL is to indicate end of args
        ... // error

    //write output
    if (sam_hdr_write(outfile, in_samhdr) < 0)
        ... // error
Refer: add_header.c

Not all type of header data can be removed but where it is possible, either a
specific header line can be removed or all of a header type can be removed. To
remove a specific line, header type, unique identifier field tag and its value
to be used. To remove all lines of a type, header type and unique identifier
field tag are to be used.

    ...

    //remove specific line
    if (sam_hdr_remove_line_id(in_samhdr, header, id, idval) < 0)
        ... // error

    //remove multiple lines of a header type
    if (sam_hdr_remove_lines(in_samhdr, header, id, NULL) < 0)
        ... // error
Refer: rem_header.c

Shows the file content after removing SQ line with SN 2.
    ./rem_header ../../samtools/test/mpileup/mpileup.1.bam SQ 2

The unique identifier for the line needs to be found to update a field, though
not all types in the header may be modifiable.  The api sam_hdr_update_line
takes the unique identifier for the header line type, its value, the field
which needs to be modified and the new value with which to modify it, followed
by a NULL.
e.g. To change LN field from 2000 to 2250 in SQ line with unique identifier SN
as 'chr1', sam_hdr_update_line( header, "SQ", "SN", "chr1", "LN", "2250",
NULL).  To change PP field from ABC to DEF in PG line with ID APP.10,
sam_hdr_update_line( header, "PG", "ID", "APP.10", "PP", "DEF", NULL).

    ...
    //update with new data
    if (sam_hdr_update_line(in_samhdr, header, id, idval, tag, val, NULL) < 0) {
        printf("Failed to update data\n");
        goto end;
    }
    ...
Refer: update_header.c

Shows new sam file with 2nd SQ line having length as 38.

    ./update_header /tmp/sample.sam.gz SQ T1 LN 38


### Update alignment data

Many of the bam data fields may be updated by setting new value to appropriate
field in bam1_core_t structure and for a few, creating a new bam1_t record would
be easier than update of existing record.

    ...
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0)
    {
    ...
            case 1:// QNAME
                ret = bam_set_qname(bamdata, val);
            break;
            case 2:// FLAG
                bamdata->core.flag = atol(val) & 0xFFFF;
            break;
            case 3:// RNAME
            case 7:// RNEXT
                if ((ret = sam_hdr_name2tid(in_samhdr, val)) < 0)
                    ... // error
                if (field == 3) {
                    //reference
                    bamdata->core.tid = ret;
                } else {
                    //mate reference
                    bamdata->core.mtid = ret;
                }
            break;
            case 4:// POS
                bamdata->core.pos = atoll(val);
            break;
            case 5:// MAPQ
                bamdata->core.qual = atoi(val) & 0x0FF;
            break;
            case 6:// CIGAR
            {
                ...
                //get cigar array and set all data in new bam record
                if ((ncigar = sam_parse_cigar(val, NULL, &cigar, &size)) < 0)
                    ... // error
                if (bam_set1(newbam, bamdata->core.l_qname, bam_get_qname(bamdata), bamdata->core.flag, bamdata->core.tid,
                 bamdata->core.pos, bamdata->core.qual, ncigar, cigar, bamdata->core.mtid, bamdata->core.mpos,
                  bamdata->core.isize, bamdata->core.l_qseq, (const char*)bam_get_seq(bamdata),
                   (const char*)bam_get_qual(bamdata), bam_get_l_aux(bamdata)) < 0)
                    ... // error

                //correct sequence data as input is expected in ascii format and not as compressed inside bam!
                memcpy(bam_get_seq(newbam), bam_get_seq(bamdata), (bamdata->core.l_qseq + 1) / 2);
                //copy the aux data
                memcpy(bam_get_aux(newbam), bam_get_aux(bamdata), bam_get_l_aux(bamdata));
            ...
            break;
            case 8:// PNEXT
                bamdata->core.mpos = atoll(val);
            break;
            case 9:// TLEN
                bamdata->core.isize = atoll(val);
            break;
            case 10:// SEQ
                ...
                for( c = 0; c < i; ++c) {
                    bam_set_seqi(bam_get_seq(bamdata), c, seq_nt16_table[(unsigned char)val[c]]);
                }
            break;
            case 11:// QUAL
                ...
                for (c = 0; c < i; ++c)
                    val[c] -= 33;               //phred score from ascii value
                memcpy(bam_get_qual(bamdata), val, i);
Refer: mod_bam.c

Shows data with RNAME modified to T2.

    ./mod_bam /tmp/sample.sam ITR1 3 T2

The auxiliary data in bam1_t structure can be modified using
bam_aux_update_float, bam_aux_update_int etc. apis. If the aux field is not
present at all, it can be appended using bam_aux_append.

    ...
    //matched to qname, update aux
    if (!(data = bam_aux_get(bamdata, tag))) {
        //tag not present append
        ... // cut: computed length and val based on tag type
        if (bam_aux_append(bamdata, tag, type, length, (const uint8_t*)val))
            ... // error
    } else {
        //update the tag with newer value
        char auxtype = bam_aux_type(data);
        switch (type) {
            case 'f':
            case 'd':
                ...
                if (bam_aux_update_float(bamdata, tag, atof(val)))
                    ... // error
            case 'C':
            case 'S':
            case 'I':
                ...
                if (bam_aux_update_int(bamdata, tag, atoll(val)))
                    ... // error
            case 'Z':
                ...
                if (bam_aux_update_str(bamdata, tag, length, val))
                    ... // error
            case 'A':
                ...
                //update the char data directly on buffer
                *(data+1) = val[0];
Refer: mod_aux.c

Shows the given record's MD tag set to Test.

    ./mod_aux samtools/test/mpileup/mpileup.1.bam ERR013140.6157908 MD Z Test

The array aux fields can be updated using bam_aux_update_array api.

    ...
    if (bam_aux_update_array(bamdata, "BA", 'I', sizeof(cnt)/sizeof(cnt[0]), cnt))
        ... // error
Refer: mod_aux_ba.c

Shows the records updated with an array of integers, containing count of ACGT
and N in that order. The bases are decoded before count for the sake of
simplicity. Refer qtask_ordered.c for a better counting where decoding is made
outside the loop.

    ./mod_aux_ba samtools/test/mpileup/mpileup.1.bam


### Create an index

Indexes help to read data faster without iterating sequentially through the
file. Indexes contain the position information about alignments and that they
can be read easily. There are different type of indices, BAI, CSI, CRAI, TBI,
FAI etc. and are usually used with iterators.

Indexing of plain/textual files are not supported, compressed SAM&FASTA/Q, BAM,
and CRAM files can be indexed. CRAM files are indexed as .crai and the others
as .bai, .csi, .fai etc. Each of these types have different internal
representations of the index information. Bai uses a fixed configuration values
where as csi has them dynamically updated based on the alignment data.

Indexes can be created either with save of alignment data or explicitly by
read of existing alignment file for alignment data (SAM/BAM/CRAM). For reference
data it has to be explicitly created (FASTA).

To create index along with alignment write, the sam_idx_init api need to be
invoked before the start of alignment data write. This api takes the output
samFile pointer, header pointer, minimum shift and index file path. For BAI
index, the min shift has to be 0.

At the end of write, sam_idx_save api need to be invoked to save the index.

    ...
    //write header
    if (sam_hdr_write(outfile, in_samhdr))
        ... // error
    // initialize indexing, before start of write
    if (sam_idx_init(outfile, in_samhdr, size, fileidx))
        ... // error
        if (sam_write1(outfile, in_samhdr, bamdata) < 0)
            ... // error
    if (sam_idx_save(outfile))
        ... // error
Refer:index_write.c

Creates mpileup.1.bam and mpileup.1.bam.bai in /tmp/.

    ./idx_on_write  ../../samtools/test/mpileup/mpileup.1.bam 0 /tmp/

To create index explicitly on an existing alignment data file, the
sam_index_build api or its alike can be used. sam_index_build takes the
alignment file path, min shift for the index and creates the index file in
same path. The output name will be based on the alignment file format and min
shift passed.

The sam_index_build2 api takes the index file path as well and gives more
control than the previous one.  The sam_index_build3 api provides an option to
configure the number of threads in index creation.

Index for reference data can be created using fai_build3 api. This creates
index file with .fai extension. If the file is bgzip-ped, a .gzi file is
created as well. It takes the path to input file and that of fai and gzi files.
When fai/gzi path are NULL, they are created along with input file.
These index files will be useful for reference data access.

    ...
    if (fai_build3(filename, NULL, NULL) == -1)
        ... // error
Refer: index_fast.c

A tabix index can be created for compressed vcf/sam/bed and other data using
tbx_index_build. It is mainly used with vcf and non-sam type files.


### Read with iterators

Index file helps to read required data without sequentially accessing the file
and are required to use iterators. The interested reference, start and end
position etc. are required to read data with iterators. With index and these
information, an iterator is created and relevant alignments can be accessed by
iterating it.

The api sam_index_load and the like does the index loading. It takes input
samFile pointer and file path. It loads the index file based on the input file
name, from the same path and with implicit index file extension - cram file
with .crai and others with .bai. The sam_index_load2 api accepts explicit path
to index file, which allows loading it from a different location and explicit
extensions. The sam_index_load3 api supports download/save of the index
locally from a remote location. These apis returns NULL on failure and index
pointer on success.

The index file path can be appended to alignment file path and used as well.
In this case the paths are expected to be separated by '##idx##'.

The sam_iter_queryi or sam_iter_querys apis may be used to create an iterator
and sam_itr_next api does the alignment data retrieval. Along with retrieval
of current data, it advances the iterator to next relevant data. The
sam_iter_queryi takes the interested positions as numeric values and
sam_iter_querys takes the interested position as a string.

With sam_iter_queryi, the reference id can be the 0 based index of reference
data, -2 for unmapped alignments, -3 to start read from beginning of file, -4
to continue from current position, -5 to return nothing. Based on the
reference id given, alignment covering the given start and end positions will
be read with sam_iter_next api.

With sam_iter_querys, the reference sequence is identified with the name and
interested positions can be described with start and end separated by '-' as
string. When sequence is identified as '.', it begins from the start of file
and when it is '*', unmapped alignments are read. Reference with <name>[:],
<name>:S, <name>:S-E, <name>:-E retrieves all data, all data covering position
S onwards, all data covering position S to E, all data covering upto position
E of reference with ID <name> respectively on read using sam_iter_next.

The index and iterator created are to be destroyed once the need is over.
sam_itr_destroy and hts_idx_destroy apis does this.

    ...
    //load index file
    if (!(idx = sam_index_load2(infile, inname, idxfile)))
        ... // error
    //create iterator
    if (!(iter = sam_itr_querys(idx, in_samhdr, region)))
        ... // error

    //read using iterator
    while ((c = sam_itr_next(infile, iter, bamdata)) >= 0)
        ... // error

    if (iter)
        sam_itr_destroy(iter);
    if (idx)
        hts_idx_destroy(idx);
    ...
Refer:index_reg_read.c

With sample.sam, region as \* will show alignments with name UNMAP2 and UNMAP3

    ./read_reg /tmp/sample.sam.gz \*

With region as \., it shows all alignments

    ./read_reg /tmp/sample.sam.gz \.

With region as T1:1-4, start 1 and end 4 it shows nothing and with T1:1-5 it
shows alignment with name ITR1.

    ./read_reg /tmp/sample.sam.gz T1:1-5

With region as T2:30-100, it shows alignment with name ITR2M which refers the
reference data T2.

    ./read_reg /tmp/sample.sam.gz T2:30-100


Multiple interested regions can be specified for read using sam_itr_regarray.
It takes index path, header, count of regions and region descriptions as array
of char array / string. This array passed need to be released by the user
itself.

    ...
    //load index file, assume it to be present in same location
    if (!(idx = sam_index_load(infile, inname)))
        ... // error
    //create iterator
    if (!(iter = sam_itr_regarray(idx, in_samhdr, regions, regcnt)))
        ... // error
    if (regions) {
        //can be freed as it is no longer required
        free(regions);
        regions = NULL;
    }

    //get required area
    while ((c = sam_itr_multi_next(infile, iter, bamdata) >= 0))
        ... // process bamdata
Refer:index_multireg_read.c

With compressed sample.sam and 2 regions from reference T1 (30 to 32) and 1
region from T2 (34 onwards), alignments with name A1, B1, A2 and ITR2M would
be shown.

    ./read_multireg /tmp/sample.sam.gz 2 T1:30-32,T2:34

To use numeric indices instead of textual regions, sam_itr_regions can be used.
It takes index file path, header, count of regions and an array of region
description (hts_reglist_t*), which has the start end positions as numerals.

The index and iterators are to be destroyed using the sam_itr_destroy and
hts_idx_destroy. The hts_reglist_t* array passed is destroyed by the library
on iterator destroy. The regions array (array of char array/string) needs to be
destroyed by the user itself.

For fasta/fastq files, the index has to be loaded using fai_load3_format which
takes the file, index file names and format. With single region specification
fai_fetch64 can be used to get bases, and fai_fetchqual64 for quality in case
of fastq data. With multiple region specification, with comma separation,
faidx_fetch_seq64 and faidx_fetch_qual64 does the job. Regions has to be parsed
using fai_parse_region in case of multiregion specifications. fai_adjust_region
is used to adjust the start-end points based on available data.

Below excerpt shows fasta/q access with single and multiregions,

    ...
    //load index
    if (!(idx = fai_load3_format(inname, NULL, NULL, FAI_CREATE, fmt)))
        ... // error

    ...
    if (!usemulti) {
        //get data from single given region
        if (!(data = fai_fetch64(idx, region, &len)))
            ... // region not found

        printf("Data: %"PRId64" %s\n", len, data);
        free((void*)data);
        //get quality for fastq type
        if (fmt == FAI_FASTQ) {
            if (!(data = fai_fetchqual64(idx, region, &len)))
                ... // region not found
        ...

    } else { // usemulti
        //parse, get each region and get data for each
        while ((remaining = fai_parse_region(idx, region, &tid, &beg, &end, HTS_PARSE_LIST))) {     //here expects regions as csv
            //parsed the region, correct end points based on actual data
            if (fai_adjust_region(idx, tid, &beg, &end) == -1)
                ... // error
            //get data for given region
            if (!(data = faidx_fetch_seq64(idx, faidx_iseq(idx, tid), beg, end, &len)))
                ... // region not found

            printf("Data: %"PRIhts_pos" %s\n", len, data);
            free((void*)data);
            data = NULL;
            //get quality data for fastq
            if (fmt == FAI_FASTQ) {
                if (!(data = faidx_fetch_qual64(idx, faidx_iseq(idx, tid), beg, end, &len)))
                    ... // error
                printf("Qual: %"PRIhts_pos" %s\n", len, data);
                free((void*)data);
            ...
            region = remaining;                                     //parse remaining region defs

    ...
    if (idx) {
        fai_destroy(idx);
    ...
Refer: read_fast_index.c


### Pileup and MPileup

Pileup shows the transposed view of the SAM alignment data, i.e. it shows the
reference positions and bases which cover that position through different reads
side by side. MPileup facilitates the piling up of multiple sam files against
each other and same reference at the same time.

Mpileup has replaced the pileup. The input expects the data to be sorted by
position.

Pileup needs to be initialized with bam_pileup_init method which takes pointer
to a method, which will be called by pileup to read data from required files,
and pointer to data which might be required for this read method to do the
read operation. It returns a pointer to the pileup iterator.

User can specify methods which need to be invoked during the load and unload
of an alignment, like constructor and destructor of objects.
Bam_plp_constructor and bam_plp_destructor methods does the setup of
these methods in the pileup iterator. During invocation of these methods, the
pointer to data passed in the initialization is passed as well. If user want
to do any custom status handling or actions during load or unload, it can be
done in these methods. Alignment specific data can be created and stored in
an argument passed to the constructor and the same will be accessible during
pileup status return. The same will be accessible during destructor as well
where any deallocation can be made.

User is expected to invoke bam_plp_auto api to get the pileup status. It
returns the pileup status or NULL on end. During this all alignments are read
one by one, using the method given in initialization for data read, until one
for a new reference is found or all alignment covering a position is read. On
such condition, the pileup status is returned and the same continuous on next
bam_plp_auto call.  The pileup status returned is an array for all positions
for which the processing is completed. Along with the result, the reference
index, position in reference data and number of alignments which covers this
position are passed. User can iterate the result array and get bases from each
alignment which covers the given reference position. The alignment specific
custom data which were created in constructor function will also be available
in the result.

The bam_plp_auto api invokes the data read method to load an alignment and the
constructor method is invoked during the load. Once the end of alignment is
passed, it is removed from the processing and destructor method is invoked,
that user could do deallocations and custom actions as in load during this
time. The custom data passed during the initialization is passed to the
constructor and destructor methods during invocation.

Once the forward and reverse strands are identified, the better of the quality
is identified and used. Both reads are required for this and hence reads are
cached until its mate is read. The maximum number of reads that can be cached
is controlled by bam_plp_set_maxcnt. Reads covering a position are cached and
as soon as mate is found, quality is adjusted and is removed from cache. Reads
above the cache limit are discarded.

Once done, the pileup iterator to be discarded by sam_plp_destroy api.

    ...
    if (!(plpiter = bam_plp_init(readdata, &conf)))
        ... // error
    //set constructor destructor callbacks
    bam_plp_constructor(plpiter, plpconstructor);
    bam_plp_destructor(plpiter, plpdestructor);

    while ((plp = bam_plp_auto(plpiter, &tid, &refpos, &n))) {
        printf("%d\t%d\t", tid+1, refpos+1);
        for (j = 0; j < n; ++j) {
            //doesnt detect succeeding insertion and deletion together here, only insertion is identified
            //deletion is detected in plp->is_del as and when pos reaches the position
            //if detection ahead is required, use bam_plp_insertion here which gives deletion length along with insertion
            if (plp[j].is_del || plp[j].is_refskip) {
                printf("*");
                continue;
            }
            //start and end are displayed in UPPER and rest on LOWER
            printf("%c", plp[j].is_head ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]) :
                            (plp[j].is_tail ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]) :
                             tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)])));
            if (plp[j].indel > 0) {
                //insertions, anyway not start or end
                printf("+%d", plp[j].indel);
                for (k = 0; k < plp[j].indel; ++k) {
                    printf("%c", tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos + k + 1)]));
                }
            }
            else if (plp[j].indel < 0) {
                printf("%d", plp[j].indel);
                for (k = 0; k < -plp[j].indel; ++k) {
                    printf("?");
                }
    ...
    if (plpiter)
        bam_plp_destroy(plpiter);
    ...
Refer:pileup.c

The read method may use a simple read or it could be an advanced read using
indices, iterators and region specifications based on the need. The constructor
method may create any custom data and store it in the pointer passed to it. The
same need to be released by use on destructor method.

MPileup works same as the pileup and supports multiple inputs against the same
reference, giving side by side view of reference and alignments from different
inputs.

MPileup needs to be initialized with bam_mpileup_init method which takes
pointer to a method, which will be called by pileup to read data from required
files, and an array of pointer to data which might be required for this read
method to do the read operation. It returns a pointer to the mpileup iterator.

User can specify methods which need to be invoked during the load and unload
of an alignment, like constructor and destructor of objects.
bam_mplp_constructor and bam_mplp_destructor methods does the setup
of these methods in the pileup iterator. During invocation of these methods,
the pointer to data passed in the initialization is passed as well. If user
want to do any custom status handling or actions during load or unload, it can
be done on these methods. Alignment specific data can be created and
stored in the custom data pointer and the same will be accessible during
return of pileup status. The same will be accessible during destructor as well
where any deallocation can be made.

User is expected to invoke bam_mplp_auto api to get the pileup status. It
returns the pileup status. During this all alignments are read one by one,
using the method given in initialization for data read, until one for a new
reference is found or all alignment covering a position is read. On such
condition, the pileup status is returned and the same continuous on next
bam_mplp_auto call.

The pileup status is returned through a parameter in the method itself, is an
array for all inputs, each containing array for positions on which the
processing is completed. Along with the result, the reference index, position
in reference data and number of alignments which covers this position are
passed. User can iterate the result array and get bases from each alignment
which covers the given reference position. The alignment specific custom data
which were created in constructor function will also be available in the
result.

Once the forward and reverse strands are identified, the better of the quality
is identified and used. Both reads are required for this and hence reads are
cached until its mate is read. The maximum number of reads that can be cached
is controlled by bam_mplp_set_maxcnt. Reads covering a position are cached and
as soon as mate is found, quality is adjusted and is removed from cache. Reads
above the cache limit are discarded.

Once done, the pileup iterator to be discarded by sam_mplp_destroy api.

    ...
    if (!(mplpiter = bam_mplp_init(argc - 1, readdata, (void**) conf)))
        ... // error
    //set constructor destructor callbacks
    bam_mplp_constructor(mplpiter, plpconstructor);
    bam_mplp_destructor(mplpiter, plpdestructor);

    while (bam_mplp64_auto(mplpiter, &tid, &refpos, depth, plp) > 0) {
        printf("%d\t%"PRIhts_pos"\t", tid+1, refpos+1);

        for (input = 0; input < argc - 1; ++input) {
            for (dpt = 0; dpt  < depth[input]; ++dpt) {
                if (plp[input][dpt].is_del || plp[input][dpt].is_refskip) {
                    printf("*");
                    continue;
                }
                //start and end are displayed in UPPER and rest on LOWER
                printf("%c", plp[input][dpt].is_head ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[input][dpt].b),
                 plp[input][dpt].qpos)]) : (plp[input]->is_tail ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[input][dpt].b),
                  plp[input][dpt].qpos)]) : tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[input][dpt].b),
                   plp[input][dpt].qpos)])));
                if (plp[input][dpt].indel > 0) {
                    //insertions, anyway not start or end
                    printf("+%d", plp[input][dpt].indel);
                    for (k = 0; k < plp[input][dpt].indel; ++k) {
                        printf("%c", tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[input][dpt].b),
                         plp[input][dpt].qpos + k + 1)]));
                    }
                }
                else if (plp[input][dpt].indel < 0) {
                    printf("%d", plp[input][dpt].indel);
                    for (k = 0; k < -plp[input][dpt].indel; ++k) {
                        printf("?");
    ...
    if (mplpiter) {
        bam_mplp_destroy(mplpiter);
    }
    ...
    if (plp) {
        free(plp);
    ...
Refer:mpileup.c

This sample takes multiple sam files and shows the pileup of data side by side.

    ./mpileup /tmp/mp.bam /tmp/mp.sam


### Base modifications

The alignment data may contain base modification information as well. This
gives the base, modifications found, orientation in which it was found and the
quality for the modification. The base modification can be identified using
hts_parse_basemod api. It stores the modification details on hts_base_mod_state
and this has to be initialized using hts_base_mod_state_alloc api.

Once the modifications are identified, they can be accessed through different
ways. bam_mods_recorded api gives the modifications identified for an alignment.
Modifications can be queried for each base position iteratively using
bam_mods_at_next_pos api. Check the returned value with buffer size to see
whether the buffer is big enough to retrieve all modifications.
Instead of querying for each position, the next modified position can be
directly retrieved directly using bam_next_basemod api. An alignment can be
queried to have a specific modification using bam_mods_query_type api. At the
end of processing, the state need to be released using hts_base_mod_state_free
api.

    ...
    if (!(ms = hts_base_mod_state_alloc()))
        ... // error
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0)
    {
        ...
        if (bam_parse_basemod(bamdata, ms))
            ... // error
        bm = bam_mods_recorded(ms, &cnt);
        for (k = 0; k < cnt; ++k) {
            printf("%c", bm[k]);
        }
        printf("\n");
        hts_base_mod mod[5] = {0};  //for ATCGN
        if (opt) {
            //option 1
            for (; i < bamdata->core.l_qseq; ++i) {
                if ((r = bam_mods_at_next_pos(bamdata, ms, mod, sizeof(mod)/sizeof(mod[0]))) <= -1) {
                    printf("Failed to get modifications\n");
                    goto end;
                }
                else if (r > (sizeof(mod) / sizeof(mod[0]))) {
                    printf("More modifications than this app can handle, update the app\n");
                    goto end;
                }
                else if (!r) {
                    //no modification at this pos
                    printf("%c", seq_nt16_str[bam_seqi(data, i)]);
                }
                //modifications
                for (j = 0; j < r; ++j) {
                    printf("%c%c%c", mod[j].canonical_base, mod[j].strand ? '-' : '+', mod[j].modified_base);
    ...
        else {
            //option 2
            while ((r = bam_next_basemod(bamdata, ms, mod, sizeof(mod)/sizeof(mod[0]), &pos)) >= 0) {
                for (; i < bamdata->core.l_qseq && i < pos; ++i) {
                    printf("%c", seq_nt16_str[bam_seqi(data, i)]);
                }
                //modifications
                for (j = 0; j < r; ++j) {
                    printf("%c%c%c", mod[j].canonical_base, mod[j].strand ? '-' : '+', mod[j].modified_base);
                }
    ...
        //check last alignment's base modification
        int strand = 0, impl = 0;
        char canonical = 0, modification[] = "mhfcgebaon";      //possible modifications
        printf("\n\nLast alignment has \n");
        for (k = 0; k < sizeof(modification) - 1; ++k) {        //avoiding NUL termination
            if (bam_mods_query_type(ms, modification[k], &strand, &impl, &canonical)) {
                printf ("No modification of %c type\n", modification[k]);
            }
            else {
                printf("%s strand has %c modified with %c, can %sassume unlisted as unmodified\n", strand ? "-/bottom/reverse" :
                "+/top/forward", canonical, modification[k], impl?"" : "not " );
            }
        }
    ...
    if (ms)
        hts_base_mod_state_free(ms);
    ...
Refer:modstate.c

The modification can be accessed in pileup mode as well. bam_mods_at_qpos gives
the modification at given pileup position. Insertion and deletion to the given
position with possible modification can be retrieved using bam_plp_insertion_mod
api.

    ...
    int plpconstructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
        //when using cd, initialize and use as it will be reused after destructor
        cd->p = hts_base_mod_state_alloc();
        //parse the bam data and gather modification data from MM tags
        return (-1 == bam_parse_basemod(b, (hts_base_mod_state*)cd->p)) ? 1 : 0;
    }

    int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
        if (cd->p) {
            hts_base_mod_state_free((hts_base_mod_state *)cd->p);
            cd->p = NULL;
        }
        return 0;
    }

    int main(int argc, char *argv[])
    {
    ...
    if (!(plpiter = bam_plp_init(readdata, &conf))) {
        ... // error
    //set constructor destructor callbacks
    bam_plp_constructor(plpiter, plpconstructor);
    bam_plp_destructor(plpiter, plpdestructor);

    while ((plp = bam_plp_auto(plpiter, &tid, &refpos, &depth))) {
        memset(&mods, 0, sizeof(mods));
        printf("%d\t%d\t", tid+1, refpos+1);

        for (j = 0; j < depth; ++j) {
            dellen = 0;
            if (plp[j].is_del || plp[j].is_refskip) {
                printf("*");
                continue;
            }
            /*invoke bam mods_mods_at_qpos before bam_plp_insertion_mod that the base modification
            is retrieved before change in pileup pos thr' plp_insertion_mod call*/
            if ((modlen = bam_mods_at_qpos(plp[j].b, plp[j].qpos, plp[j].cd.p, mods, NMODS)) == -1)
                ... // error
            //use plp_insertion/_mod to get insertion and del at the same position
            if ((inslen = bam_plp_insertion_mod(&plp[j], (hts_base_mod_state*)plp[j].cd.p, &insdata, &dellen)) == -1)
                ... // error
            //start and end are displayed in UPPER and rest on LOWER, only 1st modification considered
            //base and modification
            printf("%c%c%c", plp[j].is_head ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]) :
                (plp[j].is_tail ? toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)]) :
                    tolower(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos)])),
                modlen > 0 ? mods[0].strand ? '-' : '+' : '\0', modlen > 0 ? mods[0].modified_base : '\0');
            //insertion and deletions
            if (plp[j].indel > 0) {
                //insertion
                /*insertion data from plp_insertion_mod, note this shows the quality value as well
                which is different from base and modification above;the lower case display is not attempted either*/
                printf("+%d%s", plp[j].indel, insdata.s);
                //handle deletion if any
                if (dellen) {
                    printf("-%d", dellen);
                    for (k = 0; k < dellen; ++k) {
                        printf("?");
                ...
            else if (plp[j].indel < 0) {
                //deletion
                printf("%d", plp[j].indel);
                for (k = 0; k < -plp[j].indel; ++k) {
                    printf("?");
                }
            }
    ...
Refer:pileup_mod.c


### Read selected fields

At times the whole alignment data may not be of interest and it would be
better to read required fields alone from the alignment data. CRAM file format
supports such specific data read and HTSLib provides an option to use this.
This can improve the performance on read operation.

The hts_set_opt method does the selection of specified fields. There are flags
indicating specific fields, like SAM_FLAG, SAM_SEQ, SAM_QNAME, in alignment
data and a combination of flags for the required fields can be passed with
CRAM_OPT_REQUIRED_FIELDS to this api.

    ...
    //select required field alone, this is useful for CRAM alone
    if (hts_set_opt(infile, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG) < 0)
        ... // error

    //read header
    in_samhdr = sam_hdr_read(infile);
    ...
    //read data, check flags and update count
    while ((c = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        if (bamdata->core.flag & BAM_FREAD1)
            cntread1++;
        ...
Refer: flags_htsopt_field.c


### Thread-pool to read / write

The HTSLib api supports thread pooling for better performance. There are a few
ways in which this can be used. The pool can be made specific for a file or a
generic pool can be created and shared across multiple files. Thread pool can
also be used to execute user defined tasks. The tasks are to be added to queue,
threads in pool executes them and results can be queued back if required.

To have a thread pool specific for a file, hts_set_opt api can be used with the
file pointer, HTS_OPT_NTHREADS and the number of threads to be in the pool.
Thread pool is released on closure of file. To have a thread pool which can be
shared across different files, it needs to be initialized using hts_tpool_init
api, passing number of threads as an argument. This thread pool can be
associated with a file using hts_set_opt api. The file pointer,
HTS_OPT_THREAD_POOL and the thread pool address are to be passed as arguments to
the api. The thread pool has to be released with hts_tpool_destroy.

The samples are trivial ones to showcase the usage of api. The number of threads
to use for different tasks has to be identified based on complexity and
parallelism of the task.

Below excerpt shows file specific thread pool,

    ...
    //create file specific threads
    if (hts_set_opt(infile, HTS_OPT_NTHREADS, 1) < 0 ||     //1 thread specific for reading
    hts_set_opt(outfile1, HTS_OPT_NTHREADS, 1) < 0 ||       //1 thread specific for sam write
    hts_set_opt(outfile2, HTS_OPT_NTHREADS, 2) < 0) {       //2 thread specific for bam write
        printf("Failed to set thread options\n");
        goto end;
    }
Refer: split_thread1.c

Below excerpt shows a thread pool shared across files,

    ...
    //create a pool of 4 threads
    if (!(tpool.pool = hts_tpool_init(4)))
        ... // error
    //share the pool with all the 3 files
    if (hts_set_opt(infile, HTS_OPT_THREAD_POOL, &tpool) < 0 ||
    hts_set_opt(outfile1, HTS_OPT_THREAD_POOL, &tpool) < 0 ||
    hts_set_opt(outfile2, HTS_OPT_THREAD_POOL, &tpool) < 0) {
        ... // error

    ... // do something

    //tidy up at end
    if (tpool.pool)
        hts_tpool_destroy(tpool.pool);
    ...
Refer: split_thread2.c

Note that it is important to analyze the task in hand to decide the number of
threads to be used. As an example, if the number of threads for reading is set
to 2 and bam write to 1, keeping total number of threads the same, the
performance may decrease as bam decoding is easier than encoding.

Custom task / user defined functions can be performed on data using thread pool
and for that, the task has to be scheduled to a queue. Thread pool associated
with the queue will perform the task. There can be multiple pools and queues.
The order of execution of threads are decided based on many factors and load on
each task may vary, so the completion of the tasks may not be in the order of
their queueing. The queues can be used in two different ways, one where the
result is enqueued to queue again to be read in same order as initial queueing,
second where the resuls are not enqueued and completed possibly in a different
order than initial queueing. Explicitly created threads can also be used along
with hts thread pool usage.

hts_tpool_process_init initializes the queue / process, associates a queue with
thread pool and reserves space for given number of tasks on queue. It takes a
parameter indicating whether the result need to be enqueued for retrieval or
not. If the result is enqueued, it is retrieved in the order of scheduling of
task. Another parameter sets the maximum number of slots for tasks in queue,
usually 2 times the number of threads are used. The input and output have their
own queues and they grow as required upto the max set. hts_tpool_dispatch api
enqueues the task to the queue. The api blocks when there is no space in queue.
This behavior can be controlled with hts_tpool_dispatch2 api. The queue can be
reset using hts_tpool_process_reset api where all tasks are discarded. The api
hts_tpool_dispatch3 supports configuring cleanup routines which are to be run
when reset occurs on the queue. hts_tpool_process_flush api can ensure that
all the piled up tasks are processed, a possible case when the queueing and
processing happen at different speeds. hts_tpool_process_shutdown api stops the
processing of queue.

There are a few apis which let the user to check the status of processing. The
api hts_tpool_process_empty shows whether all the tasks are completed or not.
The api hts_tpool_process_sz gives the number of tasks, at different states of
processing. The api hts_tpool_process_len gives the number of results in output
queue waiting to be collected.

The order of execution of tasks depends on the number of threads involved and
how the threads are scheduled by operating system. When the results are enqueued
back to queue, they are read in same order of enqueueing of task and in that
case the order of execution will not be noticed. When the results are not
enqueued the results are available right away and the order of execution may be
noticeable. Based on the nature of task and the need of order maintenance, users
can select either of the queueing.

Below excerpts shows the usage of queues and threads in both cases. In the 1st,
alignments are updated with an aux tag indicating GC ratio. The order of data
has to be maintained even after update, hence the result queueing is used to
ensure same order as initial. A number of alignments are bunched together and
reuse of allocated memory is made to make it perform better. A sentinel job is
used to identify the completion of all tasks at the result collection side.
    ...
    void *thread_ordered_proc(void *args)
    {
        ...
        for ( i = 0; i < bamdata->count; ++i) {
            ...
            for (pos = 0; pos < bamdata->bamarray[i]->core.l_qseq; ++pos)
                count[bam_seqi(data,pos)]++;
            ...
            gcratio = (count[2] /*C*/ + count[4] /*G*/) / (float) (count[1] /*A*/ + count[8] /*T*/ + count[2] + count[4]);

            if (bam_aux_append(bamdata->bamarray[i], "xr", 'f', sizeof(gcratio), (const uint8_t*)&gcratio) < 0) {

    ...
    void *threadfn_orderedwrite(void *args)
    {
        ...
        //get result and write; wait if no result is in queue - until shutdown of queue
        while (tdata->result == 0 &&
            (r = hts_tpool_next_result_wait(tdata->queue)) != NULL) {
            bamdata = (data*) hts_tpool_result_data(r);
            ...
            for (i = 0; i < bamdata->count; ++i) {
                if (sam_write1(tdata->outfile, tdata->samhdr, bamdata->bamarray[i]) < 0) {
                    ... // error
            ...
            hts_tpool_delete_result(r, 0);              //release the result memory
            ...

        // Shut down the process queue.  If we stopped early due to a write failure,
        // this will signal to the other end that something has gone wrong.
        hts_tpool_process_shutdown(tdata->queue);

    ...
    int main(int argc, char *argv[])
    {
        ...
        if (!(pool = hts_tpool_init(cnt)))                  //thread pool
            ... // error
        tpool.pool = pool;      //to share the pool for file read and write as well
        //queue to use with thread pool, for task and results
        if (!(queue = hts_tpool_process_init(pool, cnt * 2, 0))) {
    ...
        //share the thread pool with i/o files
        if (hts_set_opt(infile, HTS_OPT_THREAD_POOL, &tpool) < 0 ||
            hts_set_opt(outfile, HTS_OPT_THREAD_POOL, &tpool) < 0)
            ... // error
        if (pthread_create(&thread, NULL, threadfn_orderedwrite, &twritedata))
            ... // error
        while (c >= 0) {
            if (!(bamdata = getbamstorage(chunk, &bamcache)))
                ... // error
            for (cnt = 0; cnt < bamdata->maxsize; ++cnt) {
                c = sam_read1(infile, in_samhdr, bamdata->bamarray[cnt]);
                ...
                if (hts_tpool_dispatch3(pool, queue, thread_ordered_proc, bamdata,
                                        cleanup_bamstorage, cleanup_bamstorage,
                                        0) == -1)
                    ... // error
        ...
        if (queue) {
            if (-1 == c) {
                // EOF read, send a marker to tell the threadfn_orderedwrite()
                // function to shut down.
                if (hts_tpool_dispatch(pool, queue, thread_ordered_proc,
                                    NULL) == -1) {
                    ... // error
                hts_tpool_process_shutdown(queue);

        ...
        // Wait for threadfn_orderedwrite to finish.
        if (started_thread) {
            pthread_join(thread, NULL);

        ...
        if (queue) {
            // Once threadfn_orderedwrite has stopped, the queue can be
            // cleaned up.
            hts_tpool_process_destroy(queue);
        }
    ...
Refer: qtask_ordered.c

In this 2nd, the bases are counted and GC ratio of whole file is calculated.
Order in which bases are counted is not relevant and no result queue required.
The queue is created as input only.
    ...
    void *thread_unordered_proc(void *args)
    {
        ...
        for ( i = 0; i < bamdata->count; ++i) {
            data = bam_get_seq(bamdata->bamarray[i]);
            for (pos = 0; pos < bamdata->bamarray[i]->core.l_qseq; ++pos)
                counts[bam_seqi(data, pos)]++;

        ...
        //update result and add the memory block for reuse
        pthread_mutex_lock(&bamdata->cache->lock);
        for (i = 0; i < 16; i++) {
            bamdata->bases->counts[i] += counts[i];
        }

        bamdata->next = bamdata->cache->list;
        bamdata->cache->list = bamdata;
        pthread_mutex_unlock(&bamdata->cache->lock);

    ...
    int main(int argc, char *argv[])
    {
        ...
        if (!(queue = hts_tpool_process_init(pool, cnt * 2, 1)))
            ... // error
        c = 0;
        while (c >= 0) {
            ...
            for (cnt = 0; cnt < bamdata->maxsize; ++cnt) {
                c = sam_read1(infile, in_samhdr, bamdata->bamarray[cnt]);

            ...
            if (c >= -1 ) {
                ...
                if (hts_tpool_dispatch3(pool, queue, thread_unordered_proc, bamdata,
                                        cleanup_bamstorage, cleanup_bamstorage,
                                        0) == -1)
                    ... // error
        ...
        if (-1 == c) {
            // EOF read, ensure all are processed, waits for all to finish
            if (hts_tpool_process_flush(queue) == -1) {
                fprintf(stderr, "Failed to flush queue\n");
            } else { //all done
                //refer seq_nt16_str to find position of required bases
                fprintf(stdout, "GCratio: %f\nBase counts:\n",
                    (gccount.counts[2] /*C*/ + gccount.counts[4] /*G*/) / (float)
                        (gccount.counts[1] /*A*/ + gccount.counts[8] /*T*/ +
                            gccount.counts[2] + gccount.counts[4]));
        ...
        if (queue) {
            hts_tpool_process_destroy(queue);
        }
Refer: qtask_unordered.c

## More Information

### CRAM reference files

The cram reference data is required for the read of sequence data in CRAM
format. The sequence data file may have it as embedded or as a reference to
the actual file. When it is a reference, it is downloaded locally, in the
cache directory for later usage. It will be stored in a directory structure
based on the MD5 checksum in the cache directory.

Each chromosome in a reference file gets saved as a separate file with md5sum
as its path and name. The initial 4 numerals make the directory name and rest
as the file name (<cache path>/<1st 2 of md5sum>/<2nd 2 of md5sum>/<rest of
md5sum>).

The download would be attempted from standard location, EBI ENA
(https://www.ebi.ac.uk/ena).


### Bam1_t

This structure holds the sequence data in BAM format. There are fixed and
variable size fields, basic and extended information on sequence
data. Variable size data and extended information are kept together in a
buffer, named data in bam1_t. Fields in the member named core, bam1_core_t,
and a few macros together support the storage and handling of the whole
sequence data.

- core has a link to reference as a 0 based index in field tid. The mate /
  reverse strand's link to reference is given by mtid.

- Field pos and mpos gives the position in reference to which the sequence and
  its mate / reverse strand match.

- Field flag gives the properties of the given alignment. It shows the
  alignment's orientation, mate status, read order etc.

- Field qual gives the quality of the alignment read.

- l_qname gives the length of the name of the alignment / read, l_extranul gives
  the extra space used internally in the data field.

- l_qseq gives the length of the alignment / read in the data field.

-- n_cigar gives the number of CIGAR operations for the given alignment.

- isize gives the insert size of the read / alignment.

The bases in sequence data are stored by compressing 2 bases together in a
byte.  When the reverse flag is set, the base data is reversed and
complemented from the actual read (i.e. if the forward read is ACTG, the
reverse read to be CAGT; it will be stored in SAM format with reversed and
complemented format as ACTG with reverse flag set).

Macros bam_get_qname, bam_get_seq, bam_get_qual, bam_get_aux, bam_get_l_aux,
bam_seqi etc access the data field and retrieve the required data. The aux
macros support the retrieval of auxiliary data from the data field.


### Sam_hdr_t

This structure holds the header information. This holds the number of targets
/ SQ lines in the file, each one's length, name and reference count to this
structure. It also has this information in an internal data structure for
easier access of each field of this data.

When this data is shared or assigned to another variable of a different scope
or purpose, the reference count needs to be incremented to ensure that it is
valid till the end of the variable's scope. sam_hdr_incr_ref and it needs to
be destroyed as many times with sam_hdr_destroy api.


### Index

Indices need the data to be sorted by position.  They can be of different
types with extension .bai, .csi or .tbi for compressed SAM/BAM/VCF files and
.crai for CRAM files.  The index name can be passed along with the alignment
file itself by appending a specific character sequence. The apis can detect this
sequence and extract the index path. ##idx## is the sequence which separates
the file path and index path.


### Data files

The data files can be a local file, a network file, a file accessible through
the web or in cloud storage like google and amazon.  The data files can be
represented with URIs like file://, file://localhost/.., ,ftp://..,
gs+http[s].., s3+http[s]://

