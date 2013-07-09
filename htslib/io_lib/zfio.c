#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <unistd.h>

#include "io_lib/os.h"
#include "io_lib/zfio.h"

/* ------------------------------------------------------------------------ */
/* Some wrappers around FILE * vs gzFile *, allowing for either */

/*
 * gzopen() works on both compressed and uncompressed data, but it has
 * a significant performance hit even for uncompressed data (tested as
 * 25s using FILE* to 46s via gzOpen and 66s via gzOpen when gzipped).
 *
 * Hence we use our own wrapper 'zfp' which is a FILE* when uncompressed
 * and gzFile* when compressed. This also means we could hide bzopen in
 * there too if desired.
 */

off_t zftello(zfp *zf) {
    return zf->fp ? ftello(zf->fp) : -1;
}

int zfseeko(zfp *zf, off_t offset, int whence) {
    return zf->fp ? fseeko(zf->fp, offset, whence) : -1;
}


/*
 * A wrapper for either fgets or gzgets depending on what has been
 * opened.
 */
char *zfgets(char *line, int size, zfp *zf) {
    if (zf->fp)
	return fgets(line, size, zf->fp);
    else
	return gzgets(zf->gz, line, size);
}

/*
 * A wrapper for either fputs or gzputs depending on what has been
 * opened.
 */
int zfputs(char *line, zfp *zf) {
    if (zf->fp)
	return fputs(line, zf->fp);
    else
	return gzputs(zf->gz, line) ? 0 : EOF;
}

/*
 * Peeks at and returns the next character without consuming it from the
 * input. (Ie a combination of getc and ungetc).
 */
int zfpeek(zfp *zf) {
    int c;

    if (zf->fp) {
	c = getc(zf->fp);
	if (c != EOF)
	    ungetc(c, zf->fp);
    } else {
	c = gzgetc(zf->gz);
	if (c != EOF)
	    gzungetc(c, zf->gz);
    }

    return c;
}

/* A replacement for either feof of gzeof */
int zfeof(zfp *zf) {
    return zf->fp ? feof(zf->fp) : gzeof(zf->gz);
}

/* A replacement for either fopen or gzopen */
zfp *zfopen(const char *path, const char *mode) {
    char path2[1024];
    zfp *zf;

    if (!(zf = (zfp *)malloc(sizeof(*zf))))
	return NULL;
    zf->fp = NULL;
    zf->gz = NULL;

    /* Try normal fopen */
    if (mode[0] != 'z' && mode[1] != 'z' &&
	NULL != (zf->fp = fopen(path, mode))) {
	unsigned char magic[2];
	if (2 != fread(magic, 1, 2, zf->fp))
	    return NULL;
	if (!(magic[0] == 0x1f &&
	      magic[1] == 0x8b)) {
	    fseeko(zf->fp, 0, SEEK_SET);
	    return zf;
	}

	fclose(zf->fp);
	zf->fp = NULL;
    }

#ifdef HAVE_POPEN
    /*
     * I've no idea why, by gzgets is VERY slow, maybe because it handles
     * arbitrary seeks.
     * popen to gzip -cd is 3 times faster though.
     */
    if (*mode == 'w') {
    } else {
	if (access(path, R_OK) == 0) {
	    sprintf(path2, "gzip -cd < %.*s", 1000, path);
	    if (NULL != (zf->fp = popen(path2, "r")))
		return zf;
	}
	
	sprintf(path2, "gzip -cd < %.*s.gz", 1000, path);
	if (NULL != (zf->fp = popen(path2, "r")))
	    return zf;

	printf("Failed on %s\n", path);
    } else {
	sprintf(path2, "gzip > %.*s", 1000, path);
	if (NULL != (zf->fp = popen(path2, "w")))
	    return zf;
	}
	
	printf("Failed on %s\n", path);
    }
#else
    /* Gzopen instead */
    if ((zf->gz = gzopen(path, mode)))
	return zf;

    sprintf(path2, "%.*s.gz", 1020, path);
    if ((zf->gz = gzopen(path2, mode)))
	return zf;
#endif

    perror(path);

    free(zf);
    return NULL;
}

int zfclose(zfp *zf) {
    int r = (zf->fp) ? fclose(zf->fp) : gzclose(zf->gz);
    free(zf);
    return r;
}
