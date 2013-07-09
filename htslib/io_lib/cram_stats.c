/*
 * Author: James Bonfield, Wellcome Trust Sanger Institute. 2013
 */

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <ctype.h>

#include "io_lib/cram.h"
#include "io_lib/os.h"

cram_stats *cram_stats_create(void) {
    return calloc(1, sizeof(cram_stats));
}

void cram_stats_add(cram_stats *st, int32_t val) {
    st->nsamp++;

    //assert(val >= 0);

    if (val < MAX_STAT_VAL && val >= 0) {
	st->freqs[val]++;
    } else {
	HashItem *hi;

	if (!st->h) {
	    st->h = HashTableCreate(2048, HASH_DYNAMIC_SIZE|HASH_NONVOLATILE_KEYS|HASH_INT_KEYS);
	}

	if ((hi = HashTableSearch(st->h, (char *)(size_t)val, 4))) {
	    hi->data.i++;
	} else {
	    HashData hd;
	    hd.i = 1;
	    HashTableAdd(st->h, (char *)(size_t)val, 4, hd, NULL);
	}
    }
}

void cram_stats_del(cram_stats *st, int32_t val) {
    st->nsamp--;

    //assert(val >= 0);

    if (val < MAX_STAT_VAL && val >= 0) {
	st->freqs[val]--;
	assert(st->freqs[val] >= 0);
    } else if (st->h) {
	HashItem *hi;

	if ((hi = HashTableSearch(st->h, (char *)(size_t)val, 4))) {
	    if (--hi->data.i == 0)
		HashTableDel(st->h, hi, 0);
	} else {
	    fprintf(stderr, "Failed to remove val %d from cram_stats\n", val);
	    st->nsamp++;
	}
    } else {
	fprintf(stderr, "Failed to remove val %d from cram_stats\n", val);
	st->nsamp++;
    }
}

void cram_stats_dump(cram_stats *st) {
    int i;
    fprintf(stderr, "cram_stats:\n");
    for (i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	fprintf(stderr, "\t%d\t%d\n", i, st->freqs[i]);
    }
    if (st->h) {
	HashIter *iter=  HashTableIterCreate();
	HashItem *hi;

	while ((hi = HashTableIterNext(st->h, iter))) {
	    fprintf(stderr, "\t%d\t%d\n", (int)(size_t)hi->key,
		    (int)hi->data.i);
	}
	HashTableIterDestroy(iter);
    }
}

#if 1
/* Returns the number of bits set in val; it the highest bit used */
static int nbits(int v) {
    static const int MultiplyDeBruijnBitPosition[32] = {
	1, 10, 2, 11, 14, 22, 3, 30, 12, 15, 17, 19, 23, 26, 4, 31,
	9, 13, 21, 29, 16, 18, 25, 8, 20, 28, 24, 7, 27, 6, 5, 32
    };

    v |= v >> 1; // first up to set all bits 1 after the first 1 */
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;

    // DeBruijn magic to find top bit
    return MultiplyDeBruijnBitPosition[(uint32_t)(v * 0x07C4ACDDU) >> 27];
}
#endif

/*
 * Computes entropy from integer frequencies for various encoding methods and
 * picks the best encoding.
 *
 * FIXME: we could reuse some of the code here for the actual encoding
 * parameters too. Eg the best 'k' for SUBEXP or the code lengths for huffman.
 *
 * Returns the best codec to use.
 */
enum cram_encoding cram_stats_encoding(cram_fd *fd, cram_stats *st) {
    enum cram_encoding best_encoding = E_NULL;
    int best_size = INT_MAX, bits;
    int nvals, i, ntot = 0, max_val = 0, min_val = INT_MAX, k;
    int *vals = NULL, *freqs = NULL, vals_alloc = 0, *codes;

    //cram_stats_dump(st);

    /* Count number of unique symbols */
    for (nvals = i = 0; i < MAX_STAT_VAL; i++) {
	if (!st->freqs[i])
	    continue;
	if (nvals >= vals_alloc) {
	    vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
	    vals  = realloc(vals,  vals_alloc * sizeof(int));
	    freqs = realloc(freqs, vals_alloc * sizeof(int));
	    if (!vals || !freqs)
		return E_HUFFMAN; // Cannot do much else atm
	}
	vals[nvals] = i;
	freqs[nvals] = st->freqs[i];
	ntot += freqs[nvals];
	if (max_val < i) max_val = i;
	if (min_val > i) min_val = i;
	nvals++;
    }
    if (st->h) {
	HashIter *iter=  HashTableIterCreate();
	HashItem *hi;
	int i;

	while ((hi = HashTableIterNext(st->h, iter))) {
	    if (nvals >= vals_alloc) {
		vals_alloc = vals_alloc ? vals_alloc*2 : 1024;
		vals  = realloc(vals,  vals_alloc * sizeof(int));
		freqs = realloc(freqs, vals_alloc * sizeof(int));
		if (!vals || !freqs)
		    return E_HUFFMAN; // Cannot do much else atm
	    }
	    i = (size_t)hi->key;
	    vals[nvals]=i;
	    freqs[nvals] = hi->data.i;
	    ntot += freqs[nvals];
	    if (max_val < i) max_val = i;
	    if (min_val > i) min_val = i;
	    nvals++;
	}
	HashTableIterDestroy(iter);
    }

    st->nvals = nvals;
    assert(ntot == st->nsamp);

    if (nvals <= 1) {
	free(vals);
	free(freqs);
	return E_HUFFMAN;
    }

    /*
     * Avoid complex stats for now, just do heuristic of HUFFMAN for small
     * alphabets and BETA for anything large.
     */
    free(vals); free(freqs);
    return nvals < 200 ? E_HUFFMAN : E_BETA;

    /* We only support huffman now anyway... */
    //free(vals); free(freqs); return E_HUFFMAN;

    if (fd->verbose > 1)
	fprintf(stderr, "Range = %d..%d, nvals=%d, ntot=%d\n",
		min_val, max_val, nvals, ntot);

    /* Theoretical entropy */
    {
	double dbits = 0;
	for (i = 0; i < nvals; i++) {
	    dbits += freqs[i] * log((double)freqs[i]/ntot);
	}
	dbits /= -log(2);
	if (fd->verbose > 1)
	    fprintf(stderr, "Entropy = %f\n", dbits);
    }

    /* Beta */
    bits = nbits(max_val - min_val) * ntot;
    if (fd->verbose > 1)
	fprintf(stderr, "BETA    = %d\n", bits);
    if (best_size > bits)
	best_size = bits, best_encoding = E_BETA;

#if 0
    /* Unary */
    if (min_val >= 0) {
	for (bits = i = 0; i < nvals; i++)
	    bits += freqs[i]*(vals[i]+1);
	if (fd->verbose > 1)
	    fprintf(stderr, "UNARY   = %d\n", bits);
	if (best_size > bits)
	    best_size = bits, best_encoding = E_NULL; //E_UNARY;
    }

    /* Gamma */
    for (bits = i = 0; i < nvals; i++)
	bits += ((nbits(vals[i]-min_val+1)-1) + nbits(vals[i]-min_val+1)) * freqs[i];
    if (fd->verbose > 1)
	fprintf(stderr, "GAMMA   = %d\n", bits);
    if (best_size > bits)
	best_size = bits, best_encoding = E_GAMMA;

    /* Subexponential */
    for (k = 0; k < 10; k++) {
	for (bits = i = 0; i < nvals; i++) {
	    if (vals[i]-min_val < (1<<k))
		bits += (1 + k)*freqs[i];
	    else
		bits += (nbits(vals[i]-min_val)*2-k)*freqs[i];
	}

	if (fd->verbose > 1)
	    fprintf(stderr, "SUBEXP%d = %d\n", k, bits);
	if (best_size > bits)
	    best_size = bits, best_encoding = E_SUBEXP;
    }
#endif

    /* byte array len */

    /* byte array stop */

    /* External? Guesswork! */

    /* Huffman */
//    qsort(freqs, nvals, sizeof(freqs[0]), sort_freqs);
//    for (i = 0; i < nvals; i++) {
//	fprintf(stderr, "%d = %d\n", i, freqs[i]);
//	vals[i] = 0;
//    }

    /* Grow freqs to 2*freqs, to store sums */
    /* Vals holds link data */
    freqs = realloc(freqs, 2*nvals*sizeof(*freqs));
    codes = calloc(2*nvals, sizeof(*codes));
    if (!freqs || !codes)
	return E_HUFFMAN; // Cannot do much else atm

    /* Inefficient, use pointers to form chain so we can insert and maintain
     * a sorted list? This is currently O(nvals^2) complexity.
     */
    for (;;) {
	int low1 = INT_MAX, low2 = INT_MAX;
	int ind1 = 0, ind2 = 0;
	for (i = 0; i < nvals; i++) {
	    if (freqs[i] < 0)
		continue;
	    if (low1 > freqs[i]) 
		low2 = low1, ind2 = ind1, low1 = freqs[i], ind1 = i;
	    else if (low2 > freqs[i])
		low2 = freqs[i], ind2 = i;
	}
	if (low2 == INT_MAX)
	    break;

	//fprintf(stderr, "Merge ind %d (%d), %d (%d) = %d+%d, => %d=%d\n",
	//	ind1, vals[ind1], ind2, vals[ind2], low1, low2,
	//	nvals, low1+low2);

	freqs[nvals] = low1 + low2;
	codes[ind1] = nvals;
	codes[ind2] = nvals;
	freqs[ind1] *= -1;
	freqs[ind2] *= -1;
	nvals++;
    }
    nvals = nvals/2+1;

    for (i = 0; i < nvals; i++) {
	int code_len = 0;
	for (k = codes[i]; k; k = codes[k])
	    code_len++;
	codes[i] = code_len;
	freqs[i] *= -1;
	//fprintf(stderr, "%d / %d => %d\n", vals[i], freqs[i], codes[i]);
    }

    for (bits = i = 0; i < nvals; i++) {
	bits += freqs[i] * codes[i];
    }
    if (fd->verbose > 1)
	fprintf(stderr, "HUFFMAN = %d\n", bits);
    if (best_size >= bits)
	best_size = bits, best_encoding = E_HUFFMAN;
    free(codes);

    free(vals);
    free(freqs);

    return best_encoding;
}

void cram_stats_free(cram_stats *st) {
    if (st->h)
	HashTableDestroy(st->h, 0);
    free(st);
}
