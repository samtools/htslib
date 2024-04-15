/*  sam_mods.c -- Base modification handling in SAM and BAM.

    Copyright (C) 2020-2024 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#define HTS_BUILDING_LIBRARY // Enables HTSLIB_EXPORT, see htslib/hts_defs.h
#include <config.h>
#include <assert.h>

#include "htslib/sam.h"
#include "textutils_internal.h"

// ---------------------------
// Base Modification retrieval
//
// These operate by recording state in an opaque type, allocated and freed
// via the functions below.
//
// Initially we call bam_parse_basemod to process the tags and record the
// modifications in the state structure, and then functions such as
// bam_next_basemod can iterate over this cached state.

/* Overview of API.

We start by allocating an hts_base_mod_state and parsing the MM, ML and MN
tags into it.  This has optional flags controlling how we report base
modifications in "explicit" coordinates.  See below

    hts_base_mod_state *m = hts_base_mod_state_alloc();
    bam_parse_basemod2(b, m, HTS_MOD_REPORT_UNCHECKED);
    // Or: bam_parse_basemod(b, m), which is equiv to flags==0
    //... do something ...
    hts_base_mod_state_free(m);

In the default implicit MM coordinate system, any location not
reported is implicitly assumed to contain no modification.  We only
report the places we think are likely modified.

Some tools however only look for base modifications in particular
contexts, eg CpG islands.  Here we need to distinguish between
not-looked-for and looked-for-but-didn't-find.  These calls have an
explicit coordinate system, where we only know information about the
coordinates explicitly listed and everything else is considered to be
unverified.

By default we don't get reports on the other coordinates in an
explicit MM tag, but the HTS_MOD_REPORT_UNCHECKED flag will report
them (with quality HTS_MOD_UNCHECKED) meaning we can do consensus
modification analysis with accurate counting when dealing with a
mixture of explicit and implicit records.


We have different ways of processing the base modifications.  We can
iterate either mod-by-mod or position-by-position, or we can simply
query a specific coordinate as may be done when processing a pileup.

To check for base modifications as a specific location within a
sequence we can use bam_mods_at_qpos.  This provides complete random
access within the MM string.  However currently this is inefficiently
implemented so should only be used for occasional analysis or as a way
to start iterating at a specific location.  It modifies the state
position, so after the first use we can then switch to
bam_mods_at_next_pos to iterate position by position from then on.

    hts_base_mod mods[10];
    int n = bam_mods_at_qpos(b, pos, m, mods, 10);

For base by base, we have bam_mods_at_next_pos.  This strictly starts
at the first base and reports entries one at a time.  It's more
efficient than a loop repeatedly calling ...at-pos.

    hts_base_mod mods[10];
    int n = bam_mods_at_next_pos(b, m, mods, 10);
    for (int i = 0; i < n; i++) {
        // report mod i of n
    }

Iterating over modifications instead of coordinates is simpler and
more efficient as it skips reporting of unmodified bases.  This is
done with bam_next_basemod.

    hts_base_mod mods[10];
    while ((n=bam_next_basemod(b, m, mods, 10, &pos)) > 0) {
        for (j = 0; j < n; j++) {
            // Report 'n'th mod at sequence position 'pos'
        }
    }

There are also functions that query meta-data about the MM line rather
than per-site information.

bam_mods_recorded returns an array of ints holding the +ve code ('m')
or -ve CHEBI numeric values.

    int ntypes, *types = bam_mods_recorded(m, &ntype);

We can then query a specific modification type to get further
information on the strand it is operating on, whether it has implicit
or explicit coordinates, and what it's corresponding canonical base it
is (The "C" in "C+m").  bam_mods_query_type does this by code name,
while bam_mods_queryi does this by numeric i^{th} type (from 0 to ntype-1).

    bam_mods_query_type(m, 'c', &strand, &implicit, &canonical);
    bam_mods_queryi(m, 2, &strand, &implicit, &canonical);

*/

/*
 * Base modification are stored in MM/Mm tags as <mod_list> defined as
 *
 * <mod_list>        ::= <mod_chain><mod_list> | ""
 * <mod_chain>       ::= <canonical_base><strand><mod-list><delta-list>
 *
 * <canonical_base>  ::= "A" | "C" | "G" | "T" | "N".
 *
 * <strand>          ::= "+" | "-".
 *
 * <mod-list>        ::= <simple-mod-list> | <ChEBI-code>
 * <simple-mod-list> ::= <simple-mod><simple-mod-list> | <simple-mod>
 * <ChEBI-code>      ::= <integer>
 * <simple-mod>      ::= <letter>
 *
 * <delta-list>      ::= "," <integer> <delta-list> | ";"
 *
 * We do not allocate additional memory other than the fixed size
 * state, thus we track up to 256 pointers to different locations
 * within the MM and ML tags.  Each pointer is for a distinct
 * modification code (simple or ChEBI), meaning some may point to the
 * same delta-list when multiple codes are combined together
 * (e.g. "C+mh,1,5,18,3;").  This is the MM[] array.
 *
 * Each numeric in the delta-list is tracked in MMcount[], counted
 * down until it hits zero in which case the next delta is fetched.
 *
 * ML array similarly holds the locations in the quality (ML) tag per
 * type, but these are interleaved so C+mhfc,10,15 will have 4 types
 * all pointing to the same delta position, but in ML we store
 * Q(m0)Q(h0)Q(f0)Q(c0) followed by Q(m1)Q(h1)Q(f1)Q(c1).  This ML
 * also has MLstride indicating how many positions along ML to jump
 * each time we consume a base. (4 in our above example, but usually 1
 * for the simple case).
 *
 * One complexity of the base modification system is that mods are
 * always stored in the original DNA orientation.  This is so that
 * tools that may reverse-complement a sequence (eg "samtools fastq -T
 * MM,ML") can pass through these modification tags irrespective of
 * whether they have any knowledge of their internal workings.
 *
 * Because we don't wish to allocate extra memory, we cannot simply
 * reverse the MM and ML tags.  Sadly this means we have to manage the
 * reverse complementing ourselves on-the-fly.
 * For reversed reads we start at the right end of MM and no longer
 * stop at the semicolon.  Instead we use MMend[] array to mark the
 * termination point.
 */
#define MAX_BASE_MOD 256
struct hts_base_mod_state {
    int type[MAX_BASE_MOD];     // char or minus-CHEBI
    int canonical[MAX_BASE_MOD];// canonical base, as seqi (1,2,4,8,15)
    char strand[MAX_BASE_MOD];  // strand of modification; + or -
    int MMcount[MAX_BASE_MOD];  // no. canonical bases left until next mod
    char *MM[MAX_BASE_MOD];     // next pos delta (string)
    char *MMend[MAX_BASE_MOD];  // end of pos-delta string
    uint8_t *ML[MAX_BASE_MOD];  // next qual
    int MLstride[MAX_BASE_MOD]; // bytes between quals for this type
    int implicit[MAX_BASE_MOD]; // treat unlisted positions as non-modified?
    int seq_pos;                // current position along sequence
    int nmods;                  // used array size (0 to MAX_BASE_MOD-1).
    uint32_t flags;             // Bit-field: see HTS_MOD_REPORT_UNCHECKED
};

hts_base_mod_state *hts_base_mod_state_alloc(void) {
    return calloc(1, sizeof(hts_base_mod_state));
}

void hts_base_mod_state_free(hts_base_mod_state *state) {
    free(state);
}

/*
 * Count frequency of A, C, G, T and N canonical bases in the sequence
 */
static void seq_freq(const bam1_t *b, int freq[16]) {
    int i;

    memset(freq, 0, 16*sizeof(*freq));
    uint8_t *seq = bam_get_seq(b);
    for (i = 0; i < b->core.l_qseq; i++)
        freq[bam_seqi(seq, i)]++;
    freq[15] = b->core.l_qseq; // all bases count as N for base mods
}

//0123456789ABCDEF
//=ACMGRSVTWYHKDBN  aka seq_nt16_str[]
//=TGKCYSBAWRDMHVN  comp1ement of seq_nt16_str
//084C2A6E195D3B7F
static int seqi_rc[] = { 0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15 };

/*
 * Parse the MM and ML tags to populate the base mod state.
 * This structure will have been previously allocated via
 * hts_base_mod_state_alloc, but it does not need to be repeatedly
 * freed and allocated for each new bam record. (Although obviously
 * it requires a new call to this function.)
 *
 * Flags are copied into the state and used to control reporting functions.
 * Currently the only flag is HTS_MOD_REPORT_UNCHECKED, to control whether
 * explicit "C+m?" mods report quality HTS_MOD_UNCHECKED for the bases
 * outside the explicitly reported region.
 */
int bam_parse_basemod2(const bam1_t *b, hts_base_mod_state *state,
                       uint32_t flags) {
    // Reset position, else upcoming calls may fail on
    // seq pos - length comparison
    state->seq_pos = 0;
    state->nmods = 0;
    state->flags = flags;

    // Read MM and ML tags
    uint8_t *mm = bam_aux_get(b, "MM");
    if (!mm) mm = bam_aux_get(b, "Mm");
    if (!mm)
        return 0;
    if (mm[0] != 'Z') {
        hts_log_error("%s: MM tag is not of type Z", bam_get_qname(b));
        return -1;
    }

    uint8_t *mi = bam_aux_get(b, "MN");
    if (mi && bam_aux2i(mi) != b->core.l_qseq && b->core.l_qseq) {
        // bam_aux2i with set errno = EINVAL and return 0 if the tag
        // isn't integer, but 0 will be a seq-length mismatch anyway so
        // triggers an error here too.
        hts_log_error("%s: MM/MN data length is incompatible with"
                      " SEQ length", bam_get_qname(b));
        return -1;
    }

    uint8_t *ml = bam_aux_get(b, "ML");
    if (!ml) ml = bam_aux_get(b, "Ml");
    if (ml && (ml[0] != 'B' || ml[1] != 'C')) {
        hts_log_error("%s: ML tag is not of type B,C", bam_get_qname(b));
        return -1;
    }
    uint8_t *ml_end = ml ? ml+6 + le_to_u32(ml+2) : NULL;
    if (ml) ml += 6;

    // Aggregate freqs of ACGTN if reversed, to get final-delta (later)
    int freq[16];
    if (b->core.flag & BAM_FREVERSE)
        seq_freq(b, freq);

    char *cp = (char *)mm+1;
    int mod_num = 0;
    int implicit = 1;
    while (*cp) {
        for (; *cp; cp++) {
            // cp should be [ACGTNU][+-]([a-zA-Z]+|[0-9]+)[.?]?(,\d+)*;
            unsigned char btype = *cp++;

            if (btype != 'A' && btype != 'C' &&
                btype != 'G' && btype != 'T' &&
                btype != 'U' && btype != 'N')
                return -1;
            if (btype == 'U') btype = 'T';

            btype = seq_nt16_table[btype];

            // Strand
            if (*cp != '+' && *cp != '-')
                return -1; // malformed
            char strand = *cp++;

            // List of modification types
            char *ms = cp, *me; // mod code start and end
            char *cp_end = NULL;
            int chebi = 0;
            if (isdigit_c(*cp)) {
                chebi = strtol(cp, &cp_end, 10);
                cp = cp_end;
                ms = cp-1;
            } else {
                while (*cp && isalpha_c(*cp))
                    cp++;
                if (*cp == '\0')
                    return -1;
            }

            me = cp;

            // Optional explicit vs implicit marker
            implicit = 1;
            if (*cp == '.') {
                // default is implicit = 1;
                cp++;
            } else if (*cp == '?') {
                implicit = 0;
                cp++;
            } else if (*cp != ',' && *cp != ';') {
                // parse error
                return -1;
            }

            long delta;
            int n = 0; // nth symbol in a multi-mod string
            int stride = me-ms;
            int ndelta = 0;

            if (b->core.flag & BAM_FREVERSE) {
                // We process the sequence in left to right order,
                // but delta is successive count of bases to skip
                // counting right to left.  This also means the number
                // of bases to skip at left edge is unrecorded (as it's
                // the remainder).
                //
                // To output mods in left to right, we step through the
                // MM list in reverse and need to identify the left-end
                // "remainder" delta.
                int total_seq = 0;
                for (;;) {
                    cp += (*cp == ',');
                    if (*cp == 0 || *cp == ';')
                        break;

                    delta = strtol(cp, &cp_end, 10);
                    if (cp_end == cp) {
                        hts_log_error("%s: Hit end of MM tag. Missing "
                                      "semicolon?", bam_get_qname(b));
                        return -1;
                    }

                    cp = cp_end;
                    total_seq += delta+1;
                    ndelta++;
                }
                delta = freq[seqi_rc[btype]] - total_seq; // remainder
            } else {
                delta = *cp == ','
                    ? strtol(cp+1, &cp_end, 10)
                    : 0;
                if (!cp_end) {
                    // empty list
                    delta = INT_MAX;
                    cp_end = cp;
                }
            }
            // Now delta is first in list or computed remainder,
            // and cp_end is either start or end of the MM list.
            while (ms < me) {
                state->type     [mod_num] = chebi ? -chebi : *ms;
                state->strand   [mod_num] = (strand == '-');
                state->canonical[mod_num] = btype;
                state->MLstride [mod_num] = stride;
                state->implicit [mod_num] = implicit;

                if (delta < 0) {
                    hts_log_error("%s: MM tag refers to bases beyond sequence "
                                  "length", bam_get_qname(b));
                    return -1;
                }
                state->MMcount  [mod_num] = delta;
                if (b->core.flag & BAM_FREVERSE) {
                    state->MM   [mod_num] = me+1;
                    state->MMend[mod_num] = cp_end;
                    state->ML   [mod_num] = ml ? ml+n +(ndelta-1)*stride: NULL;
                } else {
                    state->MM   [mod_num] = cp_end;
                    state->MMend[mod_num] = NULL;
                    state->ML   [mod_num] = ml ? ml+n : NULL;
                }

                if (++mod_num >= MAX_BASE_MOD) {
                    hts_log_error("%s: Too many base modification types",
                                  bam_get_qname(b));
                    return -1;
                }
                ms++; n++;
            }

            // Skip modification deltas
            if (ml) {
                if (b->core.flag & BAM_FREVERSE) {
                    ml += ndelta*stride;
                } else {
                    while (*cp && *cp != ';') {
                        if (*cp == ',')
                            ml+=stride;
                        cp++;
                    }
                }
                if (ml > ml_end) {
                    hts_log_error("%s: Insufficient number of entries in ML "
                                  "tag", bam_get_qname(b));
                    return -1;
                }
            } else {
                // cp_end already known if FREVERSE
                if (cp_end && (b->core.flag & BAM_FREVERSE))
                    cp = cp_end;
                else
                    while (*cp && *cp != ';')
                        cp++;
            }
            if (!*cp) {
                hts_log_error("%s: Hit end of MM tag. Missing semicolon?",
                              bam_get_qname(b));
                return -1;
            }
        }
    }
    if (ml && ml != ml_end) {
        hts_log_error("%s: Too many entries in ML tag", bam_get_qname(b));
        return -1;
    }

    state->nmods = mod_num;

    return 0;
}

int bam_parse_basemod(const bam1_t *b, hts_base_mod_state *state) {
    return bam_parse_basemod2(b, state, 0);
}

/*
 * Fills out mods[] with the base modifications found.
 * Returns the number found (0 if none), which may be more than
 * the size of n_mods if more were found than reported.
 * Returns <= -1 on error.
 *
 * This always marches left to right along sequence, irrespective of
 * reverse flag or modification strand.
 */
int bam_mods_at_next_pos(const bam1_t *b, hts_base_mod_state *state,
                         hts_base_mod *mods, int n_mods) {
    if (b->core.flag & BAM_FREVERSE) {
        if (state->seq_pos < 0)
            return -1;
    } else {
        if (state->seq_pos >= b->core.l_qseq)
            return -1;
    }

    int i, j, n = 0;
    unsigned char base = bam_seqi(bam_get_seq(b), state->seq_pos);
    state->seq_pos++;
    if (b->core.flag & BAM_FREVERSE)
        base = seqi_rc[base];

    for (i = 0; i < state->nmods; i++) {
        int unchecked = 0;
        if (state->canonical[i] != base && state->canonical[i] != 15/*N*/)
            continue;

        if (state->MMcount[i]-- > 0) {
            if (!state->implicit[i] &&
                (state->flags & HTS_MOD_REPORT_UNCHECKED))
                unchecked = 1;
            else
                continue;
        }

        char *MMptr = state->MM[i];
        if (n < n_mods) {
            mods[n].modified_base = state->type[i];
            mods[n].canonical_base = seq_nt16_str[state->canonical[i]];
            mods[n].strand = state->strand[i];
            mods[n].qual = unchecked
                ? HTS_MOD_UNCHECKED
                : (state->ML[i] ? *state->ML[i] : HTS_MOD_UNKNOWN);
        }
        n++;

        if (unchecked)
            continue;

        if (state->ML[i])
            state->ML[i] += (b->core.flag & BAM_FREVERSE)
                ? -state->MLstride[i]
                : +state->MLstride[i];

        if (b->core.flag & BAM_FREVERSE) {
            // process MM list backwards
            char *cp;
            if (state->MMend[i]-1 < state->MM[i]) {
                // Should be impossible to hit if coding is correct
                hts_log_error("Assert failed while processing base modification states");
                return -1;
            }
            for (cp = state->MMend[i]-1; cp != state->MM[i]; cp--)
                if (*cp == ',')
                    break;
            state->MMend[i] = cp;
            if (cp != state->MM[i])
                state->MMcount[i] = strtol(cp+1, NULL, 10);
            else
                state->MMcount[i] = INT_MAX;
        } else {
            if (*state->MM[i] == ',')
                state->MMcount[i] = strtol(state->MM[i]+1, &state->MM[i], 10);
            else
                state->MMcount[i] = INT_MAX;
        }

        // Multiple mods at the same coords.
        for (j=i+1; j < state->nmods && state->MM[j] == MMptr; j++) {
            if (n < n_mods) {
                mods[n].modified_base = state->type[j];
                mods[n].canonical_base = seq_nt16_str[state->canonical[j]];
                mods[n].strand = state->strand[j];
                mods[n].qual = state->ML[j] ? *state->ML[j] : -1;
            }
            n++;
            state->MMcount[j] = state->MMcount[i];
            state->MM[j]      = state->MM[i];
            if (state->ML[j])
                state->ML[j] += (b->core.flag & BAM_FREVERSE)
                    ? -state->MLstride[j]
                    : +state->MLstride[j];
        }
        i = j-1;
    }

    return n;
}

/*
 * Return data at the next modified location.
 *
 * bam_mods_at_next_pos does quite a bit of work, so we don't want to
 * repeatedly call it for every location until we find a mod.  Instead
 * we check how many base types we can consume before the next mod,
 * and scan through the sequence looking for them.  Once we're at that
 * site, we defer back to bam_mods_at_next_pos for the return values.
 */
int bam_next_basemod(const bam1_t *b, hts_base_mod_state *state,
                     hts_base_mod *mods, int n_mods, int *pos) {
    // Look through state->MMcount arrays to see when the next lowest is
    // per base type;
    int next[16], freq[16] = {0}, i;
    memset(next, 0x7f, 16*sizeof(*next));
    const int unchecked = state->flags & HTS_MOD_REPORT_UNCHECKED;
    if (b->core.flag & BAM_FREVERSE) {
        for (i = 0; i < state->nmods; i++) {
            if (unchecked && !state->implicit[i])
                next[seqi_rc[state->canonical[i]]] = 1;
            else if (next[seqi_rc[state->canonical[i]]] > state->MMcount[i])
                next[seqi_rc[state->canonical[i]]] = state->MMcount[i];
        }
    } else {
        for (i = 0; i < state->nmods; i++) {
            if (unchecked && !state->implicit[i])
                next[state->canonical[i]] = 0;
            else if (next[state->canonical[i]] > state->MMcount[i])
                next[state->canonical[i]] = state->MMcount[i];
       }
    }

    // Now step through the sequence counting off base types.
    for (i = state->seq_pos; i < b->core.l_qseq; i++) {
        unsigned char bc = bam_seqi(bam_get_seq(b), i);
        if (next[bc] <= freq[bc] || next[15] <= freq[15])
            break;
        freq[bc]++;
        if (bc != 15) // N
            freq[15]++;
    }
    *pos = state->seq_pos = i;

    if (b->core.flag & BAM_FREVERSE) {
        for (i = 0; i < state->nmods; i++)
            state->MMcount[i] -= freq[seqi_rc[state->canonical[i]]];
    } else {
        for (i = 0; i < state->nmods; i++)
            state->MMcount[i] -= freq[state->canonical[i]];
    }

    if (b->core.l_qseq && state->seq_pos >= b->core.l_qseq &&
        !(b->core.flag & BAM_FREVERSE)) {
        // Spots +ve orientation run-overs.
        // The -ve orientation is spotted in bam_parse_basemod2
        int i;
        for (i = 0; i < state->nmods; i++) {
            // Check if any remaining items in MM after hitting the end
            // of the sequence.
            if (state->MMcount[i] < 0x7f000000 ||
                (*state->MM[i]!=0 && *state->MM[i]!=';')) {
                hts_log_warning("MM tag refers to bases beyond sequence length");
                return -1;
            }
        }
        return 0;
    }

    int r = bam_mods_at_next_pos(b, state, mods, n_mods);
    return r > 0 ? r : 0;
}

/*
 * As per bam_mods_at_next_pos, but at a specific qpos >= the previous qpos.
 * This can only march forwards along the read, but can do so by more than
 * one base-pair.
 *
 * This makes it useful for calling from pileup iterators where qpos may
 * start part way through a read for the first occurrence of that record.
 */
int bam_mods_at_qpos(const bam1_t *b, int qpos, hts_base_mod_state *state,
                    hts_base_mod *mods, int n_mods) {
    // FIXME: for now this is inefficient in implementation.
    int r = 0;
    while (state->seq_pos <= qpos)
        if ((r = bam_mods_at_next_pos(b, state, mods, n_mods)) < 0)
            break;

    return r;
}

/*
 * Returns the list of base modification codes provided for this
 * alignment record as an array of character codes (+ve) or ChEBI numbers
 * (negative).
 *
 * Returns the array, with *ntype filled out with the size.
 *         The array returned should not be freed.
 *         It is a valid pointer until the state is freed using
 *         hts_base_mod_free().
 */
int *bam_mods_recorded(hts_base_mod_state *state, int *ntype) {
    *ntype = state->nmods;
    return state->type;
}

/*
 * Returns data about a specific modification type for the alignment record.
 * Code is either positive (eg 'm') or negative for ChEBI numbers.
 *
 * Return 0 on success or -1 if not found.  The strand, implicit and canonical
 * fields are filled out if passed in as non-NULL pointers.
 */
int bam_mods_query_type(hts_base_mod_state *state, int code,
                        int *strand, int *implicit, char *canonical) {
    // Find code entry
    int i;
    for (i = 0; i < state->nmods; i++) {
        if (state->type[i] == code)
            break;
    }
    if (i == state->nmods)
        return -1;

    // Return data
    if (strand)    *strand    = state->strand[i];
    if (implicit)  *implicit  = state->implicit[i];
    if (canonical) *canonical = "?AC?G???T??????N"[state->canonical[i]];

    return 0;
}

/*
 * Returns data about the ith modification type for the alignment record.
 *
 * Return 0 on success or -1 if not found.  The strand, implicit and canonical
 * fields are filled out if passed in as non-NULL pointers.
 */
int bam_mods_queryi(hts_base_mod_state *state, int i,
                    int *strand, int *implicit, char *canonical) {
    if (i < 0 || i >= state->nmods)
        return -1;

    // Return data
    if (strand)    *strand    = state->strand[i];
    if (implicit)  *implicit  = state->implicit[i];
    if (canonical) *canonical = "?AC?G???T??????N"[state->canonical[i]];

    return 0;
}
