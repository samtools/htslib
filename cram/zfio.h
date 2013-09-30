#ifndef _ZFIO_H_
#define _ZFIO_H_

#include <stdio.h>
#include <zlib.h>

/*
 * Either a gzFile or a FILE.
 */
typedef struct {
    FILE   *fp;
    gzFile  gz;
} zfp;

off_t zftello(zfp *zf);
int zfseeko(zfp *zf, off_t offset, int whence);
char *zfgets(char *line, int size, zfp *zf);
int zfputs(char *line, zfp *zf);
zfp *zfopen(const char *path, const char *mode);
int zfclose(zfp *zf);
int zfpeek(zfp *zf);
int zfeof(zfp *zf);

#endif /* _ZFIO_H_ */
