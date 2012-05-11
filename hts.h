#ifndef HTS_H
#define HTS_H

#include <stdint.h>

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/************
 * File I/O *
 ************/

typedef struct {
	uint32_t is_bin:1, is_write:1, is_be:1, dummy:29;
	kstring_t line;
	char *fn_aux;
	void *fp; // file pointer; actual type depending on is_bin and is_write
} htsFile;

extern int hts_verbose;

#ifdef __cplusplus
extern "C" {
#endif

	htsFile *hts_open(const char *fn, const char *mode, const char *fn_aux);
	void hts_close(htsFile *fp);

#ifdef __cplusplus
}
#endif

#endif
