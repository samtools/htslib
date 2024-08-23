/* bgzip.c -- Block compression/decompression utility.

   Copyright (C) 2008, 2009 Broad Institute / Massachusetts Institute of Technology
   Copyright (C) 2010, 2013-2019, 2021-2024 Genome Research Ltd.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notices and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdarg.h>
#include <getopt.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <sys/time.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/hfile.h"

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  include <sys/utime.h>
#endif

static const int WINDOW_SIZE = BGZF_BLOCK_SIZE;

static void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2) error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}

static int ask_yn(void)
{
    char line[1024];
    if (fgets(line, sizeof line, stdin) == NULL)
        return 0;
    return line[0] == 'Y' || line[0] == 'y';
}

static int confirm_overwrite(const char *fn)
{
    int save_errno = errno;
    int ret = 0;

    if (isatty(STDIN_FILENO)) {
        fprintf(stderr, "[bgzip] %s already exists; do you wish to overwrite (y or n)? ", fn);
        if (ask_yn()) ret = 1;
    }

    errno = save_errno;
    return ret;
}

static int known_extension(const char *ext)
{
    static const char *known[] = {
        "gz", "bgz", "bgzf",
        NULL
    };

    const char **p;
    for (p = known; *p; p++)
        if (strcasecmp(ext, *p) == 0) return 1;
    return 0;
}

static int confirm_filename(int *is_forced, const char *name, const char *ext)
{
    if (*is_forced) {
        (*is_forced)--;
        return 1;
    }

    if (!isatty(STDIN_FILENO))
        return 0;

    fprintf(stderr, "[bgzip] .%s is not a known extension; do you wish to decompress to %s (y or n)? ", ext, name);
    return ask_yn();
}

/* getfilespec - get file status data
   path        - file path for which status to be retrieved
   status      - pointer to status structure in which the data to be stored
   returns 0 on success and -1 on failure
*/
static int getfilespec(const char *path, struct stat *status)
{
    if (!path || !status) {     //invalid
        return -1;
    }
    if (!strcmp(path, "-")) {   //cant get / set for stdin/out, return success
        return 0;
    }
    if (stat(path, status) < 0) {
        return -1;
    }
    return 0;
}

/* setfilespec - set file status data
   path        - file path for which status to be set
   status      - pointer to status structure in which the data is present
   returns 0 on success and -1 on failure
   sets only the time as of now.
*/
static int setfilespec(const char *path, const struct stat *status)
{
    if (!path || !status) {     //invalid
        return -1;
    }
    if (!strcmp(path, "-")) {   //cant get / set for stdin/out, return success
        return 0;
    }

#ifdef _WIN32
    struct _utimbuf tval;
    //time upto sec - access & modification time
    tval.actime = status->st_atime;
    tval.modtime = status->st_mtime;
    if (_utime(path, &tval) < 0) {
        fprintf(stderr, "[bgzip] Failed to set file specifications.\n");
        return -1;
    }
#else
    struct timeval tval[2];
    memset(&tval[0], 0, sizeof(tval));
    //time upto sec - access time
    tval[0].tv_sec = status->st_atime;
    //time upto sec - modification time
    tval[1].tv_sec = status->st_mtime;
    if (utimes(path, &tval[0]) < 0) {
        fprintf(stderr, "[bgzip] Failed to set file specifications.\n");
        return -1;
    }
#endif //_WIN32
    return 0;
}


static int check_name_and_extension(char *name, int *forced) {
    size_t pos;
    char *ext;

    for (pos = strlen(name); pos > 0; --pos)
        if (name[pos] == '.' || name[pos] == '/') break;

    if (pos == 0 || name[pos] != '.') {
        fprintf(stderr, "[bgzip] can't find an extension in %s -- please rename\n", name);
        return 1;
    }

    name[pos] = '\0';
    ext = &name[pos+1];

    if (!(known_extension(ext) || confirm_filename(forced, name, ext))) {
        fprintf(stderr, "[bgzip] unknown extension .%s -- declining to decompress to %s\n", ext, name);
        return 2;                            //explicit N, continue and return 2
    }

    return 0;
}


static int bgzip_main_usage(FILE *fp, int status)
{
    fprintf(fp, "\n");
    fprintf(fp, "Version: %s\n", hts_version());
    fprintf(fp, "Usage:   bgzip [OPTIONS] [FILE] ...\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "   -b, --offset INT           decompress at virtual file pointer (0-based uncompressed offset)\n");
    fprintf(fp, "   -c, --stdout               write on standard output, keep original files unchanged\n");
    fprintf(fp, "   -d, --decompress           decompress\n");
    fprintf(fp, "   -f, --force                overwrite files without asking\n");
    fprintf(fp, "   -g, --rebgzip              use an index file to bgzip a file\n");
    fprintf(fp, "   -h, --help                 give this help\n");
    fprintf(fp, "   -i, --index                compress and create BGZF index\n");
    fprintf(fp, "   -I, --index-name FILE      name of BGZF index file [file.gz.gzi]\n");
    fprintf(fp, "   -k, --keep                 don't delete input files during operation\n");
    fprintf(fp, "   -l, --compress-level INT   Compression level to use when compressing; 0 to 9, or -1 for default [-1]\n");
    fprintf(fp, "   -o, --output FILE          write to file, keep original files unchanged\n");
    fprintf(fp, "   -r, --reindex              (re)index compressed file\n");
    fprintf(fp, "   -s, --size INT             decompress INT bytes (uncompressed size)\n");
    fprintf(fp, "   -t, --test                 test integrity of compressed file\n");
    fprintf(fp, "       --binary               Don't align blocks with text lines\n");
    fprintf(fp, "   -@, --threads INT          number of compression threads to use [1]\n");
    return status;
}

int main(int argc, char **argv)
{
    int c, compress, compress_level = -1, pstdout, is_forced, test, index = 0, rebgzip = 0, reindex = 0, keep, binary;
    BGZF *fp;
    char *buffer;
    long start, end, size;
    struct stat filestat;
    char *statfilename = NULL;
    char *index_fname = NULL, *write_fname = NULL;
    int threads = 1, isstdin = 0, usedstdout = 0, ret = 0, exp_out_open = 0, f_dst = -1;

    static const struct option loptions[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"offset", required_argument, NULL, 'b'},
        {"stdout", no_argument, NULL, 'c'},
        {"decompress", no_argument, NULL, 'd'},
        {"force", no_argument, NULL, 'f'},
        {"index", no_argument, NULL, 'i'},
        {"index-name", required_argument, NULL, 'I'},
        {"compress-level", required_argument, NULL, 'l'},
        {"reindex", no_argument, NULL, 'r'},
        {"rebgzip",no_argument,NULL,'g'},
        {"size", required_argument, NULL, 's'},
        {"threads", required_argument, NULL, '@'},
        {"test", no_argument, NULL, 't'},
        {"version", no_argument, NULL, 1},
        {"keep", no_argument, NULL, 'k'},
        {"binary", no_argument, NULL, 2},
        {"output", required_argument, NULL, 'o'},
        {NULL, 0, NULL, 0}
    };

    compress = 1; pstdout = 0; start = 0; size = -1; end = -1; is_forced = 0; test = 0; keep = 0; binary = 0;
    while((c  = getopt_long(argc, argv, "cdh?fb:@:s:iI:l:grtko:",loptions,NULL)) >= 0){
        switch(c){
        case 'd': compress = 0; break;
        case 'c': pstdout = 1; break;
        case 'b': start = atol(optarg); compress = 0; pstdout = 1; break;
        case 's': size = atol(optarg); pstdout = 1; break;
        case 'f': is_forced++; break;
        case 'i': index = 1; break;
        case 'I': index_fname = optarg; break;
        case 'l': compress_level = atol(optarg); break;
        case 'g': rebgzip = 1; break;
        case 'r': reindex = 1; compress = 0; break;
        case '@': threads = atoi(optarg); break;
        case 't': test = 1; compress = 0; reindex = 0; break;
        case 'k': keep = 1; break;
        case 'o': write_fname = optarg; break;
        case 1:
            printf(
"bgzip (htslib) %s\n"
"Copyright (C) 2024 Genome Research Ltd.\n", hts_version());
            return EXIT_SUCCESS;
        case  2:  binary = 1; break;
        case 'h': return bgzip_main_usage(stdout, EXIT_SUCCESS);
        case '?': return bgzip_main_usage(stderr, EXIT_FAILURE);
        }
    }
    if (size >= 0) end = start + size;
    if (end >= 0 && end < start) {
        fprintf(stderr, "[bgzip] Illegal region: [%ld, %ld]\n", start, end);
        return 1;
    }

    if ( (index || reindex) && rebgzip )
    {
        fprintf(stderr, "[bgzip] Can't produce a index and rebgzip simultaneously\n");
        return 1;
    }
    if ( rebgzip && !index_fname )
    {
        fprintf(stderr, "[bgzip] Index file name expected with rebgzip.  See -I option.\n");
        return 1;
    }
    /* avoid -I / indexfile with multiple inputs while index/reindex. these wont be set during
    read/decompress and are not considered even if set */
    if ( (index || reindex) && !write_fname && index_fname && argc - optind > 1) {
        fprintf(stderr, "[bgzip] Cannot specify index filename with multiple data file on index, reindex.\n");
        return 1;
    }

    if (write_fname) {
        if (pstdout) {
            fprintf(stderr, "[bgzip] Cannot write to %s and stdout at the same time.\n", write_fname);
            return 1;
        } else if (strncmp(write_fname, "-", strlen(write_fname)) == 0) {
            // stdout has special handling so treat as -c
            pstdout = 1;
            write_fname = NULL;
        }
    }

    do {
        isstdin = optind >= argc ? 1 : !strcmp("-", argv[optind]);          //using stdin or not?
        /* when a named output file is not used, stdout is in use when explicitly
        selected or when stdin in is in use, it needs to be closed
        explicitly to get all io errors*/

        if (!write_fname)
            usedstdout |= isstdin || pstdout || test;

        statfilename = NULL;

        if (compress == 1) {
            hFILE* f_src = NULL;
            char out_mode[3] = "w\0";
            char out_mode_exclusive[4] = "wx\0";

            if (compress_level < -1 || compress_level > 9) {
                fprintf(stderr, "[bgzip] Invalid compress-level: %d\n", compress_level);
                return 1;
            }
            if (compress_level >= 0) {
                out_mode[1] = compress_level + '0';
                out_mode_exclusive[2] = compress_level + '0';
            }
            if (!(f_src = hopen(!isstdin ? argv[optind] : "-", "r"))) {
                fprintf(stderr, "[bgzip] %s: %s\n", strerror(errno), isstdin ? "stdin" : argv[optind]);
                return 1;
            }

            if (write_fname) {
                if (!exp_out_open) {  // only open this file once for writing, close at the end
                    if ((fp = bgzf_open(write_fname, out_mode)) == NULL) {
                        fprintf(stderr, "[bgzip] can't create %s: %s\n", write_fname, strerror(errno));
                        return 1;
                    } else {
                        exp_out_open = 1;
                    }
                }
            } else if ( argc>optind && !isstdin )            //named input file that isn't an explicit "-"
            {
                if (pstdout)
                    fp = bgzf_open("-", out_mode);
                else
                {
                    char *name = malloc(strlen(argv[optind]) + 5);
                    strcpy(name, argv[optind]);
                    strcat(name, ".gz");
                    fp = bgzf_open(name, is_forced? out_mode : out_mode_exclusive);
                    if (fp == NULL && errno == EEXIST) {
                        if (confirm_overwrite(name)) {
                            fp = bgzf_open(name, out_mode);
                        }
                        else {
                            ret = 2;                        //explicit N - no overwrite, continue and return 2
                            hclose_abruptly(f_src);
                            free(name);
                            continue;
                        }
                    }
                    if (fp == NULL) {
                        fprintf(stderr, "[bgzip] can't create %s: %s\n", name, strerror(errno));
                        free(name);
                        return 1;
                    }
                    statfilename = name;
                }
            }
            else if (!pstdout && isatty(fileno((FILE *)stdout)) )
                return bgzip_main_usage(stderr, EXIT_FAILURE);
            else if ( index && !index_fname )
            {
                fprintf(stderr, "[bgzip] Index file name expected when writing to stdout\n");
                return 1;
            }
            else
                fp = bgzf_open("-", out_mode);

            if ( index ) bgzf_index_build_init(fp);
            if (threads > 1)
                bgzf_mt(fp, threads, 256);

            buffer = malloc(WINDOW_SIZE);
            if (!buffer) {
                if (statfilename) {
                    free(statfilename);
                }
                return 1;
            }
            if (rebgzip){
                if ( bgzf_index_load(fp, index_fname, NULL) < 0 ) error("Could not load index: %s.%s\n", !isstdin ? argv[optind] : index_fname, !isstdin ? "gzi" : "");

                while ((c = hread(f_src, buffer, WINDOW_SIZE)) > 0)
                    if (bgzf_block_write(fp, buffer, c) < 0) error("Could not write %d bytes: Error %d\n", c, fp->errcode);
            }
            else {
                htsFormat fmt;
                int textual = 0;
                if (!binary
                    && hts_detect_format(f_src, &fmt) == 0
                    && fmt.compression == no_compression) {
                    switch(fmt.format) {
                    case text_format:
                    case sam:
                    case vcf:
                    case bed:
                    case fasta_format:
                    case fastq_format:
                    case fai_format:
                    case fqi_format:
                        textual = 1;
                        break;
                    default: break; // silence clang warnings
                    }
                }

                if (binary || !textual) {
                    // Binary data, either detected or explicit
                    while ((c = hread(f_src, buffer, WINDOW_SIZE)) > 0)
                        if (bgzf_write(fp, buffer, c) < 0)
                            error("Could not write %d bytes: Error %d\n",
                                c, fp->errcode);
                } else {
                    /* Text mode, try a flush after a newline */
                    int in_header = 1, n = 0, long_line = 0;
                    while ((c = hread(f_src, buffer+n, WINDOW_SIZE-n)) > 0) {
                        int c2 = c+n;
                        int flush = 0;
                        if (in_header &&
                            (long_line || buffer[0] == '@' || buffer[0] == '#')) {
                            // Scan forward to find the last header line.
                            int last_start = 0;
                            n = 0;
                            while (n < c2) {
                                if (buffer[n++] != '\n')
                                    continue;

                                last_start = n;
                                if (n < c2 &&
                                    !(buffer[n] == '@' || buffer[n] == '#')) {
                                    in_header = 0;
                                    break;
                                }
                            }
                            if (!last_start) {
                                n = c2;
                                long_line = 1;
                            } else {
                                n = last_start;
                                flush = 1;
                                long_line = 0;
                            }
                        } else {
                            // Scan backwards to find the last newline.
                            n += c; // c read plus previous n overflow
                            while (--n >= 0 && ((char *)buffer)[n] != '\n')
                                ;

                            if (n >= 0) {
                                flush = 1;
                                n++;
                            } else {
                                n = c2;
                            }
                        }

                        // Pos n is either at the end of the buffer with flush==0,
                        // or the first byte after a newline and a flush point.
                        if (bgzf_write(fp, buffer, n) < 0)
                            error("Could not write %d bytes: Error %d\n",
                                n, fp->errcode);
                        if (flush)
                            if (bgzf_flush_try(fp, 65536) < 0) {// force
                                if (statfilename) {
                                    free(statfilename);
                                }
                                return -1;
                            }

                        memmove(buffer, buffer+n, c2-n);
                        n = c2-n;
                    }

                    // Trailing data.
                    if (bgzf_write(fp, buffer, n) < 0)
                        error("Could not write %d bytes: Error %d\n",
                            n, fp->errcode);
                }
            }
            if ( index && !write_fname )
            {
                if (index_fname) {
                    if (bgzf_index_dump(fp, index_fname, NULL) < 0)
                        error("Could not write index to '%s'\n", index_fname);
                } else if (!isstdin) {
                    if (bgzf_index_dump(fp, argv[optind], ".gz.gzi") < 0)
                        error("Could not write index to '%s.gz.gzi'\n", argv[optind]);
                }
                else {
                    //stdin, cant create index file as name is not present "-.gz.gzi" not a valid one!
                    error("Can not write index for stdin data without index filename, use -I option to set index file.\n");
                }
            }

            if (!write_fname) {
                if (bgzf_close(fp) < 0)
                    error("Output close failed: Error %d\n", fp->errcode);
            }

            if (hclose(f_src) < 0)
                error("Input close failed\n");

            if (statfilename) {
                //get input file timestamp
                if (!getfilespec(argv[optind], &filestat)) {
                    //set output file timestamp
                    if (setfilespec(statfilename, &filestat) < 0) {
                        fprintf(stderr, "[bgzip] Failed to set file specification.\n");
                    }
                }
                else {
                    fprintf(stderr, "[bgzip] Failed to get file specification.\n");
                }
                free(statfilename);
            }

            if (argc > optind && !pstdout && !keep && !isstdin && !write_fname) unlink(argv[optind]);

            free(buffer);
        }
        else if ( reindex )
        {
            if ( argc>optind && !isstdin )
            {
                fp = bgzf_open(argv[optind], "r");
                if ( !fp ) error("[bgzip] Could not open file: %s\n", argv[optind]);
            }
            else
            {
                if ( !index_fname ) error("[bgzip] Index file name expected when reading from stdin\n");
                fp = bgzf_open("-", "r");
                if ( !fp ) error("[bgzip] Could not read from stdin: %s\n", strerror(errno));
            }

            buffer = malloc(BGZF_BLOCK_SIZE);
            bgzf_index_build_init(fp);
            int ret;
            while ( (ret=bgzf_read(fp, buffer, BGZF_BLOCK_SIZE))>0 ) ;
            free(buffer);
            if ( ret<0 ) error("Is the file gzipped or bgzipped? The latter is required for indexing.\n");

            if ( index_fname ) {
                if (bgzf_index_dump(fp, index_fname, NULL) < 0)
                    error("Could not write index to '%s'\n", index_fname);
            } else if (!isstdin) {
                if (bgzf_index_dump(fp, argv[optind], ".gzi") < 0)
                    error("Could not write index to '%s.gzi'\n", argv[optind]);
            }
            else {
                //stdin, cant create index file as name is not present "-.gzi" not a valid one!
                error("Can not write index for stdin data without index filename, use -I option to set index file.\n");
            }

            if ( bgzf_close(fp)<0 ) error("Close failed: Error %d\n",fp->errcode);
        }
        else
        {
            int is_forced_tmp = is_forced;

            if ( argc>optind && !isstdin )
            {
                fp = bgzf_open(argv[optind], "r");
                if (fp == NULL) {
                    fprintf(stderr, "[bgzip] Could not open %s: %s\n", argv[optind], strerror(errno));
                    return 1;
                }
                if (bgzf_compression(fp) == no_compression) {
                    fprintf(stderr, "[bgzip] %s: not a compressed file -- ignored\n", argv[optind]);
                    bgzf_close(fp);
                    return 1;
                }

                if (pstdout || test) {
                    f_dst = fileno(stdout);
                } else {
                    const int wrflags = O_WRONLY | O_CREAT | O_TRUNC;
                    char *name;
                    int check;

                    if (!(name = strdup(argv[optind]))) {
                        fprintf(stderr, "[bgzip] unable to allocate memory for output file name.\n");
                        bgzf_close(fp);
                        return 1;
                    }

                    if ((check = check_name_and_extension(name, &is_forced_tmp))) {
                        bgzf_close(fp);

                        if (check == 1) {
                            return 1;
                        } else {
                            ret = 2;
                            continue;
                        }
                    }

                    if (!exp_out_open) {
                        if (write_fname) { // only open file once and don't care about overwriting
                            is_forced_tmp = 1;
                            exp_out_open = 1;
                        }

                        f_dst = open(write_fname ? write_fname : name, is_forced_tmp? wrflags : wrflags|O_EXCL, 0666);

                        if (f_dst < 0 && errno == EEXIST) {
                            if (confirm_overwrite(name)) {
                                f_dst = open(name, wrflags, 0666);
                            }
                            else {
                                ret = 2;                        //explicit N - no overwrite, continue and return 2
                                bgzf_close(fp);
                                free(name);
                                continue;
                            }
                        }
                        if (f_dst < 0) {
                            fprintf(stderr, "[bgzip] can't create %s: %s\n", name, strerror(errno));
                            free(name);
                            return 1;
                        }
                    }

                    statfilename = name;
                }
            }
            else if (!pstdout && isatty(fileno((FILE *)stdin)) )
                return bgzip_main_usage(stderr, EXIT_FAILURE);
            else
            {
                f_dst = fileno(stdout);
                fp = bgzf_open("-", "r");
                if (fp == NULL) {
                    fprintf(stderr, "[bgzip] Could not read from stdin: %s\n", strerror(errno));
                    return 1;
                }
                if (bgzf_compression(fp) == no_compression) {
                    fprintf(stderr, "[bgzip] stdin is not compressed -- ignored\n");
                    bgzf_close(fp);
                    return 1;
                }

                if (!write_fname) {
                    f_dst = fileno(stdout);
                } else {
                    if (!exp_out_open) {
                        exp_out_open = 1;

                        f_dst = open(write_fname, O_WRONLY | O_CREAT | O_TRUNC, 0666);

                        if (f_dst < 0) {
                            fprintf(stderr, "[bgzip] can't create %s: %s\n", write_fname, strerror(errno));
                            return 1;
                        }
                    }
                }
            }

            buffer = malloc(WINDOW_SIZE);
            if ( start>0 )
            {
                if (index_fname) {
                    if ( bgzf_index_load(fp, index_fname, NULL) < 0 )
                        error("Could not load index: %s\n", index_fname);
                } else {
                    if (optind >= argc || isstdin) {
                        error("The -b option requires -I when reading from stdin "
                            "(and stdin must be seekable)\n");
                    }
                    if ( bgzf_index_load(fp, argv[optind], ".gzi") < 0 )
                        error("Could not load index: %s.gzi\n", argv[optind]);
                }
                if ( bgzf_useek(fp, start, SEEK_SET) < 0 ) error("Could not seek to %ld-th (uncompressd) byte\n", start);
            }

            if (threads > 1)
                bgzf_mt(fp, threads, 256);

    #ifdef _WIN32
            _setmode(f_dst, O_BINARY);
    #endif
            long start_reg = start, end_reg = end;
            while (1) {
                if (end < 0) c = bgzf_read(fp, buffer, WINDOW_SIZE);
                else c = bgzf_read(fp, buffer, (end - start > WINDOW_SIZE)? WINDOW_SIZE:(end - start));
                if (c == 0) break;
                if (c < 0) error("Error %d in block starting at offset %" PRId64 "(%" PRIX64 ")\n", fp->errcode, fp->block_address, fp->block_address);
                start += c;
                if ( !test && write(f_dst, buffer, c) != c ) {
    #ifdef _WIN32
                    if (GetLastError() != ERROR_NO_DATA)
    #endif
                    error("Could not write %d bytes\n", c);
                }
                if (end >= 0 && start >= end) break;
            }
            start = start_reg;
            end = end_reg;
            free(buffer);
            if (bgzf_close(fp) < 0) error("Close failed: Error %d\n",fp->errcode);

            if (statfilename) {
                if (!write_fname) {
                    //get input file timestamp
                    if (!getfilespec(argv[optind], &filestat)) {
                        //set output file timestamp
                        if (setfilespec(statfilename, &filestat) < 0) {
                            fprintf(stderr, "[bgzip] Failed to set file specification.\n");
                        }
                    }
                    else {
                        fprintf(stderr, "[bgzip] Failed to get file specification.\n");
                    }
                }

                free(statfilename);
            }

            if (argc > optind && !pstdout && !test && !keep && !isstdin && !write_fname) unlink(argv[optind]);
            if (!isstdin && !pstdout && !test && !write_fname) {
                close(f_dst);                               //close output file when it is not stdout
            }
        }
    } while (++optind < argc);

    if (usedstdout && !reindex) {
        //stdout in use, have to close explicitly to get any pending write errors
        if (fclose(stdout) != 0 && errno != EBADF) {
            fprintf(stderr, "[bgzip] Failed to close stdout, errno %d", errno);
            ret = 1;
        }
    } else if (write_fname) {
        if (compress == 1) { // close explicit output file (this is for compression)
            if (index) {
                if (index_fname) {
                    if (bgzf_index_dump(fp, index_fname, NULL) < 0)
                        error("Could not write index to '%s'\n", index_fname);
                } else {
                    if (bgzf_index_dump(fp, write_fname, ".gzi") < 0)
                        error("Could not write index to '%s.gzi'\n", write_fname);
                }
            }

            if (bgzf_close(fp) < 0)
                error("Output close failed: Error %d\n", fp->errcode);
        } else {
            close(f_dst);
        }
    }


    return ret;
}
