#ifndef _STRING_ALLOC_H_
#define _STRING_ALLOC_H_

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * A pooled string allocator intended to cut down on the
 * memory overhead of many small string allocations.
 *
 * Andrew Whitwham, September 2010.
 */

typedef struct {
    char *str;
    size_t used;
} string_t;

typedef struct {
    size_t max_length;
    size_t nstrings;
    string_t *strings;
} string_alloc_t;

string_alloc_t *string_pool_create(size_t max_length);
void string_pool_destroy(string_alloc_t *a_str);
char *string_alloc(string_alloc_t *a_str, size_t length);
char *string_dup(string_alloc_t *a_str, char *instr);
char *string_ndup(string_alloc_t *a_str, char *instr, size_t len);

#endif

#ifdef __cplusplus
}
#endif

