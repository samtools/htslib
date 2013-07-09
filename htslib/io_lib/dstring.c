#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <io_lib/dstring.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <errno.h>

#include "io_lib/dstring.h"
#include "io_lib/vlen.h"

/*
 * Allocates a new dstring, initialising it to a default str (or NULL).
 *
 * Returns dstring_t pointer on success.
 *         NULL on failure.
 */
dstring_t *dstring_create(const char *str) {
    dstring_t *ds = (dstring_t *)malloc(sizeof(*ds));
    if (!ds)
	return NULL;

    ds->str = NULL;
    ds->allocated = 0;
    ds->length = 0;

    if (str) {
	if (dstring_insert(ds, 0, str) == -1) {
	    dstring_destroy(ds);
	    return NULL;
	}
    }

    return ds;
}

/*
 * As per dstring_create(), but using str,len as the internal data.
 * Ie the caller is giving this data to the dstring object. str should
 * be a malloced pointer.
 *
 * Returns dstring_t pointer on success.
 *         NULL on failure.
 */
dstring_t *dstring_create_with(char *str, size_t len) {
    dstring_t *ds = (dstring_t *)malloc(sizeof(*ds));
    if (!ds)
	return NULL;

    ds->str = str;
    ds->allocated = ds->length = len;

    return ds;
}

/* Deallocates a dstring */
void dstring_destroy(dstring_t *ds) {
    if (ds) {
	if (ds->str)
	    free(ds->str);
	free(ds);
    }
}

/*
 * Returns a C string from a dstring. If the dstring is empty this may be
 * NULL.
 */
char *dstring_str(const dstring_t *ds) {
    return ds->str;
}

/*
 * Empties a dstring without freeing the contents. Ie sets it to contain
 * a blank string.
 */
void dstring_empty(dstring_t *ds) {
    ds->length = 0;
    if (ds->str)
	*ds->str = 0;
}

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
int dstring_resize(dstring_t *ds, size_t length) {
    char *str;

    if (length+1 <= ds->allocated)
	return 0;

    /*
     * Allocate with additional overhead so as to reduce calling this
     * to often. Increase to next power of 2.
     */
    length = pow(2, ceil(log(length+1)/log(2)));
    /* length++;*/
    str = realloc(ds->str, length);
    if (!str)
	return -1;
    else {
	ds->allocated = length;
	/* If this is first alloc, make sure we null terminate */
	if (!ds->str) {
	    str[0] = 0;
	}
	ds->str = str;
    }

    return 0;
}

/*
 * Uses vsprintf style formatting to produce a string which is then inserted
 * at a specific offset.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_vinsertf(dstring_t *ds,
		     size_t offset,
		     const char *fmt,
		     va_list args) {
    size_t est_length;
    char buf[8192], *bufp = buf;

    /*
     * Work out the expanded length, if it's small enough use a temporary
     * buffer to store the sprintf output.
     */
    est_length = vflen((char *)fmt, args);
    
    if (est_length+1 > 8192) {
	if (NULL == (bufp = (char *)malloc(est_length+1))) {
	    goto error;
	}
    }
    
    /* Produce the C string, and add it to the dstring */
    vsprintf(bufp, fmt, args);

    if (-1 == dstring_insert(ds, offset, bufp))
	goto error;

    if (bufp != buf)
	free(bufp);

    va_end(args);

    return 0;

 error:
    if (bufp && bufp != buf)
	free(bufp);

    va_end(args);

    return -1;
}

/*
 * Uses sprintf style formatting to produce a string which is then inserted
 * at a specific offset.
 *
 * Returns 0 for success
 *        -1 for failure
 */
__PRINTF_FORMAT__(3,4)
int dstring_insertf(dstring_t *ds,
		    size_t offset,
		    const char *fmt,
		    ...) {
    va_list args;
    va_start(args, fmt);

    return dstring_vinsertf(ds, offset, fmt, args);
}

/*
 * Inserts a string of a specified length to a dstring at a given offset.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_ninsert(dstring_t *ds,
		    size_t offset,
		    const char *str,
		    size_t len) {
    if (0 != DSTRING_RESIZE(ds, ds->length + len))
	return -1;

    memmove(&ds->str[offset+len], &ds->str[offset], ds->length + 1 - offset);
    memmove(&ds->str[offset], str, len);

    ds->length += len;

    return 0;
}

/*
 * Inserts a C string into a dstring at a given offset.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_insert(dstring_t *ds, size_t offset, const char *str) {
    return dstring_ninsert(ds, offset, str, strlen(str));
}

/*
 * Adds C string to the start.
 * Equivalent a dstring_insert at position zero.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_prepend(dstring_t *ds, const char *str) {
    return dstring_insert(ds, 0, str);
}

/*
 * Adds string of a specified length to the start.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_nprepend(dstring_t *ds, const char *str, size_t len) {
    return dstring_ninsert(ds, 0, str, len);
}

/*
 * Uses sprintf style formatting to produce a string which is then inserted
 * to the start of our dstring.
 *
 * Returns 0 for success
 *        -1 for failure
 */
__PRINTF_FORMAT__(2,3)
int dstring_prependf(dstring_t *ds, const char *fmt, ...) {
    
    va_list args;
    va_start(args, fmt);

    return dstring_vinsertf(ds, 0, fmt, args);
}

/*
 * Adds C string to the end.
 * Equivalent a dstring_insert at position <dstring length>.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_append(dstring_t *ds, const char *str) {
    return dstring_insert(ds, ds->length, str);
}

/*
 * Adds a C string to the end, but URL encoding it with percent rules
 * If specified 'meta' is a list of meta-characters need escaping, otherwise
 * we use standard html ones ("<>&");
 */
int dstring_append_hex_encoded(dstring_t *ds, const char *str,
			       const char *meta) {
    unsigned char escape[256];
    const unsigned char *ustr = (const unsigned char *)str;
    int i, j;
    char hex[3];

    for (i = 0; i < 256; i++) {
	if (isprint(i))
	    escape[i] = 0;
	else
	    escape[i] = 1;
    }
    escape['%'] = 1;
    if (meta) {
	for (i = 0; meta[i]; i++)
	    escape[(uc)meta[i]] = 1;
    } else {
	for (i = 0; "<>&"[i]; i++)
	    escape[(uc)"<>&"[i]] = 1;
    }

    j = 0;
    hex[0] = '%';
    do {
	i = j;
	while (ustr[i] && !escape[ustr[i]])
	    i++;

	if (i-j) {
	    if (0 != dstring_ninsert(ds, ds->length, &str[j], i-j))
		return -1;
	}

	while (ustr[i] && escape[ustr[i]]) {
	    hex[1] = "0123456789ABCDEF"[ustr[i] >>  4];
	    hex[2] = "0123456789ABCDEF"[ustr[i] & 0xf];
	    if (0 != dstring_ninsert(ds, ds->length, hex, 3))
		return -1;
	    i++;
	}
	j = i;
    } while (ustr[i]);

    return 0;
}

/*
 * Adds a single character to the end of the string.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_append_char(dstring_t *ds, char c) {
    return dstring_ninsert(ds, ds->length, &c, 1);
}

/*
 * Adds an integer to the end of the string, turning it to printable ascii.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_append_int(dstring_t *ds, int i) {
    char buf[50], *cp = buf;
    int j, k = 0;

    if (i == 0) {
	*cp++ = '0';
    } else {
	if (i < 0) {
	    *cp++ = '-';
	    i = -i;
	}

	if (i < 1000)
	    goto b1;
	if (i < 100000)
	    goto b2;
	if (i < 100000000)
	    goto b3;

	j = i / 1000000000;
	if (j || k) *cp++ = j + '0', k=1, i %= 1000000000;

	j = i / 100000000;
	if (j || k) *cp++ = j + '0', k=1, i %= 100000000;
    
    b3:
	j = i / 10000000;
	if (j || k) *cp++ = j + '0', k=1, i %= 10000000;
    
	j = i / 1000000;
	if (j || k) *cp++ = j + '0', k=1, i %= 1000000;
    
	j = i / 100000;
	if (j || k) *cp++ = j + '0', k=1, i %= 100000;
    
    b2:
	j = i / 10000;
	if (j || k) *cp++ = j + '0', k=1, i %= 10000;

	j = i / 1000;
	if (j || k) *cp++ = j + '0', k=1, i %= 1000;

    b1:
	j = i / 100;
	if (j || k) *cp++ = j + '0', k=1, i %= 100;

	j = i / 10;
	if (j || k) *cp++ = j + '0', k=1, i %= 10;

	if (i || k) *cp++ = i + '0';
    }

    return dstring_ninsert(ds, ds->length, buf, cp-buf);
}

/*
 * Adds string of a specified length to the end.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_nappend(dstring_t *ds, const char *str, size_t len) {
    if (0 != DSTRING_RESIZE(ds, ds->length + len))
	return -1;

    memcpy(&ds->str[ds->length], str, len);
    ds->length += len;

    return 0;
}

/*
 * Uses sprintf style formatting to produce a string which is then appended
 * to our dstring.
 *
 * Returns 0 for success
 *        -1 for failure
 */
__PRINTF_FORMAT__(2,3)
int dstring_appendf(dstring_t *ds, const char *fmt, ...) {
    
    va_list args;
    va_start(args, fmt);

    return dstring_vinsertf(ds, ds->length, fmt, args);
}

/*
 * Refreshes the cached dstring length.
 * Use this if you obtain a copy of the internal C string and manipulate it
 * in some way.
 */
void dstring_refresh_length(dstring_t *ds) {
    ds->length = strlen(ds->str);
}

/*
 * Returns the length of the dstring (excluding nul; like strlen).
 */
size_t dstring_length(dstring_t *ds) {
    return ds->length;
}

/*
 * Inserts a dstring into a dstring
 *
 * Returns 0 for success
 *        -1 for failure
 */
int dstring_dinsert(dstring_t *ds_to,
		    size_t offset,
		    const dstring_t *ds_from) {
    if (!ds_from || !ds_to)
	return -1;

    return dstring_insert(ds_to, offset, ds_from->str);
}

/*
 * Deletes a section from a dstring, starting at 'offset' and extending
 * for 'length' characters.
 */
void dstring_delete(dstring_t *ds, size_t offset, size_t length) {
    memmove(&ds->str[offset], &ds->str[offset+length],
	    ds->length + 1 - (offset + length));
    ds->length -= length;
}

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
		    const char *rep_str) {
    size_t rep_len = strlen(rep_str);
    
    /* Ensure our string is large enough */
    if (rep_len > length) {
	if (0 != DSTRING_RESIZE(ds, ds->length + rep_len - length))
	    return -1;
    }
    
    /* Do the replace */
    if (rep_len != length) {
	memmove(&ds->str[offset+rep_len], &ds->str[offset+length],
		ds->length + 1 - (offset + length));
    }
    memmove(&ds->str[offset], rep_str, rep_len);
    ds->length += rep_len - length;

    return 0;
}

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
		     const dstring_t *rep_with) {
    return dstring_replace(ds, offset, length, rep_with->str);
}


/*
 * Searches for the first occurance of 'search' in a dstring starting
 * at position offset (including looking at that position).
 *
 * Returns the new offset if found
 *        -1 if not.
 */
int dstring_find(dstring_t *ds,
		 size_t offset,
		 const char *search) {
    size_t i;
    size_t search_len = strlen(search);

    /* Noddy algorithm to start with; use Boyer-Moore or something if needed */
    for (i = offset; i <= ds->length; i++) {
	if (strncmp(&ds->str[i], search, search_len) == 0)
	    return i;
    }

    return -1;
}

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
			 const char *rep_with) {
    int pos;
    size_t search_len = strlen(search);

    /* Find */
    if (-1 == (pos = dstring_find(ds, offset, search)))
	return -1;

    /* And replace */
    if (0 != dstring_replace(ds, pos, search_len, rep_with))
	return -1;

    return pos;
}

/*
 * Look for 'search' starting at a specific offset. If found replace it with
 * replace. Repeat until all occurances have been replaced.
 *
 * Returns 0 for success
 *        -1 on error
 */
int dstring_find_replace_all(dstring_t *ds,
			     const char *search,
			     const char *rep_with) {
    /*
     * Most efficient mechanism is to search and replace to a separate
     * buffer so that we don't have problems with search being contained
     * within rep_with and do not have efficiency problems due to
     * excessive use of string shifting.
     */
    dstring_t *new_ds = dstring_create(NULL);
    int found_pos, current_pos = 0;
    size_t search_len = strlen(search);
    dstring_t tmp;

    if (!new_ds)
	goto error;

    while (-1 != (found_pos = dstring_find(ds, current_pos, search))) {
	if (-1 == dstring_nappend(new_ds, &ds->str[current_pos],
				  found_pos - current_pos))
	    goto error;
	if (-1 == dstring_append(new_ds, rep_with))
	    goto error;
	current_pos = found_pos + search_len;
    }
    if (-1 == dstring_append(new_ds, &ds->str[current_pos]))
	goto error;

    /*
     * Swap data structures over and free the original dstring. We can't
     * swap pointers as the calling function still needs 'ds' to be the
     * same address.
     */
    tmp = *ds;
    *ds = *new_ds;
    *new_ds = tmp;
    dstring_destroy(new_ds);

    return 0;

 error:
    if (new_ds)
	dstring_destroy(new_ds);
    return -1;
}

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
int dstring_escape_html(dstring_t *ds) {
    if (-1 == dstring_find_replace_all(ds, "&", "&amp;"))
	return -1;
    if (-1 == dstring_find_replace_all(ds, "<", "&lt;"))
	return -1;
    if (-1 == dstring_find_replace_all(ds, ">", "&gt;"))
	return -1;
    if (-1 == dstring_find_replace_all(ds, "\"", "&quot;"))
	return -1;

    return 0;
}

/*
 * Searches for URLs in text strings and converts then to html href links.
 * At present we just look for http://, https://, ftp://, file:// and
 * mailto://
 *
 * Returns 0 for success
 *        -1 on error
 */
int dstring_htmlise_links(dstring_t *ds) {
    char *links[] = {"http://", "https://", "ftp://", "file://", "mailto://"};
    size_t nlinks = sizeof(links)/sizeof(*links);
    size_t i;

    for (i = 0; i < nlinks; i++) {
	int pos = 0;
	while (-1 != (pos = dstring_find(ds, pos, links[i]))) {
	    size_t end, url_len;
	    char *str = dstring_str(ds);
	    dstring_t *url;

	    /* Got the left end; extend to find right end */
	    for (end = pos+1; str[end] && !isspace(str[end]); end++)
		;

	    /* Create a new href string */
	    if (NULL == (url = dstring_create(NULL)))
		return -1;

	    if (-1 == dstring_insertf(url, 0, "<a href=\"%.*s\">%.*s</a>",
				      (int)(end-pos), &str[pos],
				      (int)(end-pos), &str[pos])) {
		dstring_destroy(url);
		return -1;
	    }

	    url_len = dstring_length(url);
	    if (-1 == dstring_dreplace(ds, pos, end-pos, url)) {
		dstring_destroy(url);
		return -1;
	    }
	    dstring_destroy(url); 
	    pos += url_len;
	}
    }
    
    return 0;
}

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
int dstring_to_html(dstring_t *ds) {
    if (-1 == dstring_escape_html(ds))
	return -1;

    return dstring_htmlise_links(ds);
}

/*
 * Appends an system error much like perror(), using errno.
 * 'str' is added in the form "str: error_message".
 *
 * All Return 0 for success
 *           -1 for failure
 */
int dstring_perror(dstring_t *ds, const char *str) {
    return dstring_appendf(ds, "%s: %s\n", str, strerror(errno));
}

#ifdef TEST

int main(void) {
    dstring_t *ds1 = dstring_create("foo");
    dstring_t *ds2 = dstring_create(NULL);
    dstring_t *ds3 = dstring_create("blah");

    printf("ds1='%s'\n", dstring_str(ds1));
    printf("ds3='%s'\n", dstring_str(ds3));

    dstring_resize(ds2, 10);

    printf("ret=%d, ", dstring_insert(ds1, 3, "abc"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_insert(ds1, 0, "xyz"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_dinsert(ds1, 3, ds3));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_replace(ds1, 3, 4, "fish"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_replace(ds1, 3, 4, "X"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_replace(ds1, 3, 1, "YZ"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_replace(ds1, 3, 2, ""));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_replace(ds1, 3, 0, "blah"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d\n", dstring_find(ds1, 0, "ab"));
    printf("ret=%d\n", dstring_find(ds1, 0, "a"));
    printf("ret=%d\n", dstring_find(ds1, 6, "a"));
    printf("ret=%d\n", dstring_find(ds1, 11, "a"));

    printf("ret=%d, ", dstring_replace(ds1, 3, 0, "blah"));
    printf("ds1='%s'\n", dstring_str(ds1));
    printf("len=%ld, strlen=%ld\n",
	   (long)dstring_length(ds1),
	   (long)strlen(dstring_str(ds1)));

    printf("ret=%d, ", dstring_dreplace(ds1, 0, 10, ds3));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_find_replace(ds1, 0, "h", "abc"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_find_replace_all(ds1, "ab", "X"));
    printf("ds1='%s'\n", dstring_str(ds1));
    printf("ret=%d, ", dstring_find_replace_all(ds1, "b", "X"));
    printf("ret=%d, ", dstring_find_replace_all(ds1, "l", "XX"));
    printf("ret=%d, ", dstring_find_replace_all(ds1, "oo", "X"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_find_replace_all(ds1, "XX", "XXX"));
    printf("ds1='%s'\n", dstring_str(ds1));

    printf("ret=%d, ", dstring_find_replace_all(ds1, "X", ""));
    printf("ds1='%s'\n", dstring_str(ds1));

    {
	int i;
	for (i = 0; i < 8; i++) {
	    dstring_find_replace_all(ds1, "h", "hh");
	    dstring_find_replace_all(ds1, "hh", "xhhhx");
	    printf("i=%d, len=%ld\n", i, (long)strlen(dstring_str(ds1)));
	}
    }

    dstring_destroy(ds1);
    dstring_destroy(ds2);
    dstring_destroy(ds3);

    ds1 = dstring_create("xyz<http://www.mrc-lmb.cam.ac.uk/ foo&bar>=\"fish\"");
    printf("ds1='%s'\n", dstring_str(ds1));
    dstring_to_html(ds1);
    printf("ds1='%s'\n", dstring_str(ds1));

    return 0;
}

#endif
