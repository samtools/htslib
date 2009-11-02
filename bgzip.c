/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/select.h>
#include "bgzf.h"

static const int WINDOW_SIZE = 64 * 1024;

static int is_ready(int fd)
{
	fd_set fdset;
	struct timeval timeout;
	FD_ZERO(&fdset);
	FD_SET(fd, &fdset);
	timeout.tv_sec = 0; timeout.tv_usec = 1;
	return select(1, &fdset, NULL, NULL, &timeout) == 1 ? 1 : 0;
}

static int bgzip_main_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   bgzip [options] [file] ...\n\n");
	fprintf(stderr, "Options: -c      write on standard output, keep original files unchanged\n");
	fprintf(stderr, "         -d      decompress\n");
	fprintf(stderr, "         -b INT  decompress at virtual file pointer INT\n");
	fprintf(stderr, "         -s INT  decompress INT bytes in the uncompressed file\n");
	fprintf(stderr, "         -h      give this help\n");
	fprintf(stderr, "\n");
	return 1;
}

static int write_open(const char *fn, int is_forced)
{
	int fd = -1;
	char c;
	if (!is_forced) {
		if ((fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC | O_EXCL, 0666)) < 0 && errno == EEXIST) {
			fprintf(stderr, "[bgzip] %s already exists; do you wish to overwrite (y or n)? ", fn);
			scanf("%c", &c);
			if (c != 'Y' && c != 'y') {
				fprintf(stderr, "[bgzip] not overwritten\n");
				exit(1);
			}
		}
	}
	if (fd < 0) {
		if ((fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC, 0666)) < 0) {
			fprintf(stderr, "[bgzip] %s: Fail to write\n", fn);
			exit(1);
		}
	}
	return fd;
}

static void fail(BGZF* fp)
{
    fprintf(stderr, "Error: %s\n", fp->error);
    exit(1);
}

int main(int argc, char **argv)
{
	int c, compress, pstdout, is_forced;
	BGZF *fp;
	void *buffer;
	long start, end, size;

	compress = 1; pstdout = 0; start = 0; size = -1; end = -1; is_forced = 0;
	while((c  = getopt(argc, argv, "cdhfb:s:")) >= 0){
		switch(c){
		case 'h': return bgzip_main_usage();
		case 'd': compress = 0; break;
		case 'c': pstdout = 1; break;
		case 'b': start = atol(optarg); break;
		case 's': size = atol(optarg); break;
		case 'f': is_forced = 1; break;
		}
	}
	if (size >= 0) end = start + size;
	if (end >= 0 && end < start) {
		fprintf(stderr, "[bgzip] Illegal region: [%ld, %ld]\n", start, end);
		return 1;
	}
	if (compress == 1) {
		int f_src, f_dst = -1;
		if (is_ready(fileno(stdin))) pstdout = 1;
		if (argc > optind && !pstdout) {
			if ((f_src = open(argv[optind], O_RDONLY)) < 0) {
				fprintf(stderr, "[bgzip] Cannot open file: %s\n", argv[optind]);
				return 1;
			}
			if (pstdout) {
				f_dst = fileno(stdout);
			} else {
				char *name = malloc(sizeof(strlen(argv[optind]) + 5));
				strcpy(name, argv[optind]);
				strcat(name, ".gz");
				f_dst = write_open(name, is_forced);
				if (f_dst < 0) return 1;
				free(name);
			}
		} else if (pstdout) { 
			f_src = fileno(stdin);
			f_dst = fileno(stdout);
		} else return bgzip_main_usage();
		fp = bgzf_fdopen(f_dst, "w");
		buffer = malloc(WINDOW_SIZE);
		while ((c = read(f_src, buffer, WINDOW_SIZE)) > 0)
			if (bgzf_write(fp, buffer, c) < 0) fail(fp);
		// f_dst will be closed here
		if (bgzf_close(fp) < 0) fail(fp);
		if (argc > optind) unlink(argv[optind]);
		free(buffer);
		close(f_src);
		return 0;
	} else {
		int f_dst, is_stdin = 0;
		if (argc == optind) pstdout = 1;
		if (is_ready(fileno(stdin))) is_stdin = 1;
		if (argc <= optind && !is_stdin) return bgzip_main_usage();
		if (argc > optind && !pstdout) {
			char *name;
			if (strstr(argv[optind], ".gz") - argv[optind] != strlen(argv[optind]) - 3) {
				fprintf(stderr, "[bgzip] %s: unknown suffix -- ignored\n", argv[optind]);
				return 1;
			}
			name = strdup(argv[optind]);
			name[strlen(name) - 3] = '\0';
			f_dst = write_open(name, is_forced);
			free(name);
		} else f_dst = fileno(stdout);
		fp = (argc == optind)? bgzf_fdopen(fileno(stdin), "r") : bgzf_open(argv[optind], "r");
		if (fp == NULL) {
			fprintf(stderr, "[bgzip] Could not open file: %s\n", argv[optind]);
			return 1;
		}
		buffer = malloc(WINDOW_SIZE);
		if (bgzf_seek(fp, start, SEEK_SET) < 0) fail(fp);
		while (1) {
			if (end < 0) c = bgzf_read(fp, buffer, WINDOW_SIZE);
			else c = bgzf_read(fp, buffer, (end - start > WINDOW_SIZE)? WINDOW_SIZE:(end - start));
			if (c == 0) break;
			if (c < 0) fail(fp);
			start += c;
			write(f_dst, buffer, c);
			if (end >= 0 && start >= end) break;
		}
		free(buffer);
		if (bgzf_close(fp) < 0) fail(fp);
		if (!pstdout) unlink(argv[optind]);
		return 0;
	}
}
