#ifndef _DSTRING_H
#define _DSTRING_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdarg.h>

#include "io_lib/misc.h"

/*
 * Implements a simple dynamic string object.
 * Like C, offsets start from 0.
 */

typedef struct {
    char *str;		/* String ptr itself */
    size_t allocated;	/* Amount of memory malloced (including the nul) */
    size_t length;	/* Amount of memory used (excluding the nul) */
} dstring_t;


#define DSTRING_STR(ds) ((ds)->str)
#define DSTRING_LEN(ds) ((ds)->length)

/*
 * Allocates a new dstring, initialising it to a default str (or NULL).
 *
 * Returns dstring_t pointer on success.
 *         NULL on failure.
 */
dstring_t *dstring_create(const char *str);

/*
 * As per dstring_create(), but using str,len as the internal data.
 * Ie the caller is giving this data to the dstring object. str should
 * be a malloced pointer.
 *
 * Returns dstring_t pointer on success.
 *         NULL on failure.
 */
dstring_t *dstring_create_with(char *str, size_t len);

/* Deallocates a dstring */
void dstring_destroy(dstring_t *ds);


/*
 * Returns a C string from a dstring. If the dstring is empty this may be
 * NULL.
 */
char *dstring_str(const dstring_t *ds);

/*
 * Force the memory allocated for a dstring to be at least length characters
 * long. (The allocated length will include 1 more to allow for the nul
 * termination.)
 * It's possible to shrink a string too, although shrinking a string will not
 * guarantee if remains nul terminated.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_resize(dstring_t *ds, size_t length);

#define DSTRING_RESIZE(ds, len) ((ds)->allocated > (len) ? 0 : dstring_resize((ds),(len)))

/*
 * Refreshes the cached dstring length.
 * Use this if you obtain a copy of the internal C string and manipulate it
 * in some way.
 */
void dstring_refresh_length(dstring_t *ds);


/*
 * Returns the length of the dstring (excluding nul; like strlen).
 */
size_t dstring_length(dstring_t *ds);


/*
 * Insertion functions.
 * dstring_ninsert, nappend and nprepend take a string and a length (much
 * like strncmp, strncpy, etc).
 * dstring_insert, append and prepend just take a normal C string.
 * dstring_dinsert inserts one dstring into another.
 *
 * All Return 0 for success
 *           -1 for failure
 */
int dstring_insert(dstring_t *ds, size_t offset, const char *str);
int dstring_ninsert(dstring_t *ds,
		    size_t offset,
		    const char *str,
		    size_t len);
int dstring_dinsert(dstring_t *ds_to,
		    size_t offset,
		    const dstring_t *ds_from);
int dstring_vinsertf(dstring_t *ds,
		     size_t offset,
		     const char *fmt,
		     va_list args);
int dstring_insertf(dstring_t *ds, size_t offset, const char *fmt, ...) __PRINTF_FORMAT__(3,4);
int dstring_prepend(dstring_t *ds, const char *str);
int dstring_nprepend(dstring_t *ds, const char *str, size_t len);
int dstring_prependf(dstring_t *ds, const char *fmt, ...) __PRINTF_FORMAT__(2,3);
int dstring_append(dstring_t *ds, const char *str);
int dstring_nappend(dstring_t *ds, const char *str, size_t len);
int dstring_appendf(dstring_t *ds, const char *fmt, ...) __PRINTF_FORMAT__(2,3);
int dstring_append_char(dstring_t *ds, char c);
int dstring_append_int(dstring_t *ds, int i);
int dstring_append_hex_encoded(dstring_t *ds, const char *str,
			       const char *meta);

void dstring_empty(dstring_t *ds);


/*
 * Deletes a section from a dstring, starting at 'offset' and extending
 * for 'length' characters.
 */
void dstring_delete(dstring_t *ds, size_t offset, size_t length);

/*
 * Replaces a section from a dstring (at offset for length bytes) with a
 * new (C) string.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_replace(dstring_t *ds,
		    size_t offset,
		    size_t length,
		    const char *rep_str);

/*
 * Replaces a section from a dstring (at offset for length bytes) with a
 * new dstring.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_dreplace(dstring_t *ds,
		     size_t offset,
		     size_t length,
		     const dstring_t *rep_with);

/*
 * Searches for the first occurance of 'search' in a dstring starting
 * at position offset (including looking at that position).
 *
 * Returns the new offset if found
 *        -1 if not.
 */
int dstring_find(dstring_t *ds,
		 size_t offset,
		 const char *search);

/*
 * A combination of dstring_find and dstring_replace.
 * Look for 'search' starting at a specific offset. If found replace it with
 * replace.
 *
 * Returns position of replaced string if found
 *        -1 if not found or on error.
 */
int dstring_find_replace(dstring_t *ds,
			 size_t offset,
			 const char *search,
			 const char *rep_with);

/*
 * Look for 'search' starting at a specific offset. If found replace it with
 * replace. Repeat until all occurances have been replaced.
 *
 * Returns 0 for success
 *        -1 on error
 */
int dstring_find_replace_all(dstring_t *ds,
			     const char *search,
			     const char *rep_with);

/*
 * Converts a text string into a HTML version representing the same string.
 * This includes escaping any HTML meta characters and searching for URLs
 * within the string and replacing it with an HTML link (keeping the link as
 * the anchor name).
 * This is simply a wrapper joining dstring_escape_html and
 * dstring_htmlise_links.
 *
 * Returns 0 for success
 *        -1 on error
 */
int dstring_to_html(dstring_t *ds);

/*
 * Escapes HTML meta characters by replacing them with appropriate HTML
 * codes.
 * We deal with the following:
 *
 * &	&amp;
 * <	&lt;
 * >	&gt;
 * "	&quot;
 *
 * Returns 0 for success
 *        -1 on error
 */
int dstring_escape_html(dstring_t *ds);

/*
 * Searches for URLs in text strings and converts then to html href links.
 * At present we just look for http://, https://, ftp://, file:// and
 * mailto://
 *
 * Returns 0 for success
 *        -1 on error
 */
int dstring_htmlise_links(dstring_t *ds);

/*
 * Appends an system error much like perror().
 * 'str' is added in the form "str: error_message".
 *
 * All Return 0 for success
 *           -1 for failure
 */
int dstring_perror(dstring_t *ds, const char *str);

#ifdef __cplusplus
}
#endif

#endif
