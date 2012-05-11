#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bgzf.h"
#include "endian.h"
#include "hts.h"

#include "kseq.h"
KSTREAM_INIT2(, gzFile, gzread, 16384)

int hts_verbose = 3;

htsFile *hts_open(const char *fn, const char *mode, const char *fn_aux)
{
	htsFile *fp;
	fp = (htsFile*)calloc(1, sizeof(htsFile));
	fp->is_be = ed_is_big();
	if (strchr(mode, 'w')) fp->is_write = 1;
	if (strchr(mode, 'b')) fp->is_bin = 1;
	if (fp->is_bin) {
		if (fp->is_write) fp->fp = strcmp(fn, "-")? bgzf_open(fn, mode) : bgzf_dopen(fileno(stdout), mode);
		else fp->fp = strcmp(fn, "-")? bgzf_open(fn, "r") : bgzf_dopen(fileno(stdin), "r");
	} else {
		if (!fp->is_write) {
			gzFile gzfp;
			gzfp = strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(fileno(stdin), "rb");
			if (gzfp) fp->fp = ks_init(gzfp);
			if (fn_aux) fp->fn_aux = strdup(fn_aux);
		} else fp->fp = strcmp(fn, "-")? fopen(fn, "rb") : stdout;
	}
	if (fp->fp == 0) {
		if (hts_verbose >= 2)
			fprintf(stderr, "[E::%s] fail to open file '%s'\n", __func__, fn);
		free(fp->fn_aux); free(fp);
		return 0;
	}
	return fp;
}

void hts_close(htsFile *fp)
{
	if (!fp->is_bin) {
		free(fp->line.s);
		if (!fp->is_write) {
			gzFile gzfp = ((kstream_t*)fp->fp)->f;
			ks_destroy((kstream_t*)fp->fp);
			gzclose(gzfp);
			free(fp->fn_aux);
		} else fclose((FILE*)fp->fp);
	} else bgzf_close((BGZF*)fp->fp);
	free(fp);
}
