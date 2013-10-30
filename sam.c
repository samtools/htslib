#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <zlib.h>
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "cram/cram.h"
#include "hfile.h"

#include "htslib/khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int64_t)

typedef khash_t(s2i) sdict_t;

/**********************
 *** BAM header I/O ***
 **********************/

bam_hdr_t *bam_hdr_init()
{
	return (bam_hdr_t*)calloc(1, sizeof(bam_hdr_t));
}

void bam_hdr_destroy(bam_hdr_t *h)
{
	int32_t i;
	if (h == NULL) return;
	if (h->target_name) {
		for (i = 0; i < h->n_targets; ++i)
			free(h->target_name[i]);
		free(h->target_name);
		free(h->target_len);
	}
	free(h->text); free(h->cigar_tab);
	if (h->sdict) kh_destroy(s2i, (sdict_t*)h->sdict);
	free(h);
}

static bam_hdr_t *hdr_from_dict(sdict_t *d)
{
	bam_hdr_t *h;
	khint_t k;
	h = bam_hdr_init();
	h->sdict = d;
	h->n_targets = kh_size(d);
	h->target_len = (uint32_t*)malloc(4 * h->n_targets);
	h->target_name = (char**)malloc(sizeof(void*) * h->n_targets);
	for (k = kh_begin(d); k != kh_end(d); ++k) {
		if (!kh_exist(d, k)) continue;
		h->target_name[kh_val(d, k)>>32] = (char*)kh_key(d, k);
		h->target_len[kh_val(d, k)>>32]  = kh_val(d, k)<<32>>32;
		kh_val(d, k) >>= 32;
	}
	return h;
}

bam_hdr_t *bam_hdr_read(BGZF *fp)
{
	bam_hdr_t *h;
	char buf[4];
	int magic_len, has_EOF;
	int32_t i = 1, name_len;
	// check EOF
	has_EOF = bgzf_check_EOF(fp);
	if (has_EOF < 0) {
		perror("[W::sam_hdr_read] bgzf_check_EOF");
	} else if (has_EOF == 0 && hts_verbose >= 2)
		fprintf(stderr, "[W::%s] EOF marker is absent. The input is probably truncated.\n", __func__);
	// read "BAM1"
	magic_len = bgzf_read(fp, buf, 4);
	if (magic_len != 4 || strncmp(buf, "BAM\1", 4)) {
		if (hts_verbose >= 1) fprintf(stderr, "[E::%s] invalid BAM binary header\n", __func__);
		return 0;
	}
	h = bam_hdr_init();
	// read plain text and the number of reference sequences
	bgzf_read(fp, &h->l_text, 4);
	if (fp->is_be) ed_swap_4p(&h->l_text);
	h->text = (char*)malloc(h->l_text + 1);
	h->text[h->l_text] = 0; // make sure it is NULL terminated
	bgzf_read(fp, h->text, h->l_text);
	bgzf_read(fp, &h->n_targets, 4);
	if (fp->is_be) ed_swap_4p(&h->n_targets);
	// read reference sequence names and lengths
	h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
	h->target_len = (uint32_t*)calloc(h->n_targets, 4);
	for (i = 0; i != h->n_targets; ++i) {
		bgzf_read(fp, &name_len, 4);
		if (fp->is_be) ed_swap_4p(&name_len);
		h->target_name[i] = (char*)calloc(name_len, 1);
		bgzf_read(fp, h->target_name[i], name_len);
		bgzf_read(fp, &h->target_len[i], 4);
		if (fp->is_be) ed_swap_4p(&h->target_len[i]);
	}
	return h;
}

int bam_hdr_write(BGZF *fp, const bam_hdr_t *h)
{
	char buf[4];
	int32_t i, name_len, x;
	// write "BAM1"
	strncpy(buf, "BAM\1", 4);
	bgzf_write(fp, buf, 4);
	// write plain text and the number of reference sequences
	if (fp->is_be) {
		x = ed_swap_4(h->l_text);
		bgzf_write(fp, &x, 4);
		if (h->l_text) bgzf_write(fp, h->text, h->l_text);
		x = ed_swap_4(h->n_targets);
		bgzf_write(fp, &x, 4);
	} else {
		bgzf_write(fp, &h->l_text, 4);
		if (h->l_text) bgzf_write(fp, h->text, h->l_text);
		bgzf_write(fp, &h->n_targets, 4);
	}
	// write sequence names and lengths
	for (i = 0; i != h->n_targets; ++i) {
		char *p = h->target_name[i];
		name_len = strlen(p) + 1;
		if (fp->is_be) {
			x = ed_swap_4(name_len);
			bgzf_write(fp, &x, 4);
		} else bgzf_write(fp, &name_len, 4);
		bgzf_write(fp, p, name_len);
		if (fp->is_be) {
			x = ed_swap_4(h->target_len[i]);
			bgzf_write(fp, &x, 4);
		} else bgzf_write(fp, &h->target_len[i], 4);
	}
	bgzf_flush(fp);
	return 0;
}

int bam_name2id(bam_hdr_t *h, const char *ref)
{
	sdict_t *d = (sdict_t*)h->sdict;
	khint_t k;
	if (h->sdict == 0) {
		int i, absent;
		d = kh_init(s2i);
		for (i = 0; i < h->n_targets; ++i) {
			k = kh_put(s2i, d, h->target_name[i], &absent);
			kh_val(d, k) = i;
		}
		h->sdict = d;
	}
	k = kh_get(s2i, d, ref);
	return k == kh_end(d)? -1 : kh_val(d, k);
}

/*************************
 *** BAM alignment I/O ***
 *************************/

bam1_t *bam_init1()
{
	return (bam1_t*)calloc(1, sizeof(bam1_t));
}

void bam_destroy1(bam1_t *b)
{
	if (b == 0) return;
	free(b->data); free(b);
}

bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
{
	uint8_t *data = bdst->data;
	int m_data = bdst->m_data;   // backup data and m_data
	if (m_data < bsrc->l_data) { // double the capacity
		m_data = bsrc->l_data; kroundup32(m_data);
		data = (uint8_t*)realloc(data, m_data);
	}
	memcpy(data, bsrc->data, bsrc->l_data); // copy var-len data
	*bdst = *bsrc; // copy the rest
	// restore the backup
	bdst->m_data = m_data;
	bdst->data = data;
	return bdst;
}

int bam_cigar2qlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; ++k)
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			l += bam_cigar_oplen(cigar[k]);
	return l;
}

int bam_cigar2rlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; ++k)
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			l += bam_cigar_oplen(cigar[k]);
	return l;
}

int32_t bam_endpos(const bam1_t *b)
{
	if (!(b->core.flag & BAM_FUNMAP) && b->core.n_cigar > 0)
		return b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
	else
		return b->core.pos + 1;
}

static inline int aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f') return 4;
	else return 0;
}

static void swap_data(const bam1_core_t *c, int l_data, uint8_t *data)
{
	uint8_t *s;
	uint32_t *cigar = (uint32_t*)(data + c->l_qname);
	int i;
	s = data + c->n_cigar*4 + c->l_qname + c->l_qseq + (c->l_qseq + 1)/2;
	for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
	while (s < data + l_data) {
		uint8_t type;
		s += 2; // skip key
		type = toupper(*s); ++s; // skip type
		if (type == 'C' || type == 'A') ++s;
		else if (type == 'S') { ed_swap_2p(s); s += 2; }
		else if (type == 'I' || type == 'F') { ed_swap_4p(s); s += 4; }
		else if (type == 'D') { ed_swap_8p(s); s += 8; }
		else if (type == 'Z' || type == 'H') { while (*s) ++s; ++s; }
		else if (type == 'B') {
			int32_t n, Bsize = aux_type2size(*s);
			memcpy(&n, s + 1, 4);
			if (Bsize == 2) {
				for (i = 0; i < n; i += 2)
					ed_swap_2p(s + 5 + i);
			} else if (Bsize == 4) {
				for (i = 0; i < n; i += 4)
					ed_swap_4p(s + 5 + i);
			}
			ed_swap_4p(s+1); 
		}
	}
}

int bam_read1(BGZF *fp, bam1_t *b)
{
	bam1_core_t *c = &b->core;
	int32_t block_len, ret, i;
	uint32_t x[8];
	if ((ret = bgzf_read(fp, &block_len, 4)) != 4) {
		if (ret == 0) return -1; // normal end-of-file
		else return -2; // truncated
	}
	if (bgzf_read(fp, x, 32) != 32) return -3;
	if (fp->is_be) {
		ed_swap_4p(&block_len);
		for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
	}
	c->tid = x[0]; c->pos = x[1];
	c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
	c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
	c->l_qseq = x[4];
	c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];
	b->l_data = block_len - 32;
	if (b->m_data < b->l_data) {
		b->m_data = b->l_data;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	if (bgzf_read(fp, b->data, b->l_data) != b->l_data) return -4;
	//b->l_aux = b->l_data - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
	if (fp->is_be) swap_data(c, b->l_data, b->data);
	return 4 + block_len;
}

int bam_write1(BGZF *fp, const bam1_t *b)
{
	const bam1_core_t *c = &b->core;
	uint32_t x[8], block_len = b->l_data + 32, y;
	int i;
	x[0] = c->tid;
	x[1] = c->pos;
	x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | c->l_qname;
	x[3] = (uint32_t)c->flag<<16 | c->n_cigar;
	x[4] = c->l_qseq;
	x[5] = c->mtid;
	x[6] = c->mpos;
	x[7] = c->isize;
	bgzf_flush_try(fp, 4 + block_len);
	if (fp->is_be) {
		for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
		y = block_len;
		bgzf_write(fp, ed_swap_4p(&y), 4);
		swap_data(c, b->l_data, b->data);
	} else bgzf_write(fp, &block_len, 4);
	bgzf_write(fp, x, 32);
	bgzf_write(fp, b->data, b->l_data);
	if (fp->is_be) swap_data(c, b->l_data, b->data);
	return 4 + block_len;
}

/********************
 *** BAM indexing ***
 ********************/

static hts_idx_t *bam_index(BGZF *fp, int min_shift)
{
	int n_lvls, i, fmt;
	bam1_t *b;
	hts_idx_t *idx;
	bam_hdr_t *h;
	h = bam_hdr_read(fp);
	if (min_shift > 0) {
		int64_t max_len = 0, s;
		for (i = 0; i < h->n_targets; ++i)
			if (max_len < h->target_len[i]) max_len = h->target_len[i];
		max_len += 256;
		for (n_lvls = 0, s = 1<<min_shift; max_len > s; ++n_lvls, s <<= 3);
		fmt = HTS_FMT_CSI;
	} else min_shift = 14, n_lvls = 5, fmt = HTS_FMT_BAI;
	idx = hts_idx_init(h->n_targets, fmt, bgzf_tell(fp), min_shift, n_lvls);
	bam_hdr_destroy(h);
	b = bam_init1();
	while (bam_read1(fp, b) >= 0) {
		int l, ret;
		l = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
		if (l == 0) l = 1; // no zero-length records
		ret = hts_idx_push(idx, b->core.tid, b->core.pos, b->core.pos + l, bgzf_tell(fp), !(b->core.flag&BAM_FUNMAP));
		if (ret < 0) break; // unsorted
	}
	hts_idx_finish(idx, bgzf_tell(fp));
	bam_destroy1(b);
	return idx;
}

int bam_index_build(const char *fn, int min_shift)
{
	hts_idx_t *idx;
	BGZF *fp;
	if ((fp = bgzf_open(fn, "r")) == 0) return -1;
	idx = bam_index(fp, min_shift);
	bgzf_close(fp);
	hts_idx_save(idx, fn, min_shift > 0? HTS_FMT_CSI : HTS_FMT_BAI);
	hts_idx_destroy(idx);
	return 0;
}

int bam_readrec(BGZF *fp, void *null, bam1_t *b, int *tid, int *beg, int *end)
{
	int ret;
	if ((ret = bam_read1(fp, b)) >= 0) {
		*tid = b->core.tid; *beg = b->core.pos;
		*end = b->core.pos + (b->core.n_cigar? bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) : 1);
	}
	return ret;
}

/**********************
 *** SAM header I/O ***
 **********************/

#include "htslib/kseq.h"
#include "htslib/kstring.h"

bam_hdr_t *sam_hdr_parse(int l_text, const char *text)
{
	const char *q, *r, *p;
	khash_t(s2i) *d;
	d = kh_init(s2i);
	for (p = text; *p; ++p) {
		if (strncmp(p, "@SQ", 3) == 0) {
			char *sn = 0;
			int ln = -1;
			for (q = p + 4;; ++q) {
				if (strncmp(q, "SN:", 3) == 0) {
					q += 3;
					for (r = q; *r != '\t' && *r != '\n'; ++r);
					sn = (char*)calloc(r - q + 1, 1);
					strncpy(sn, q, r - q);
					q = r;
				} else if (strncmp(q, "LN:", 3) == 0)
					ln = strtol(q + 3, (char**)&q, 10);
				while (*q != '\t' && *q != '\n') ++q;
				if (*q == '\n') break;
			}
			p = q;
			if (sn && ln >= 0) {
				khint_t k;
				int absent;
				k = kh_put(s2i, d, sn, &absent);
				if (!absent) {
					if (hts_verbose >= 2)
						fprintf(stderr, "[W::%s] duplicated sequence '%s'\n", __func__, sn);
					free(sn);
				} else kh_val(d, k) = (int64_t)(kh_size(d) - 1)<<32 | ln;
			}
		}
		while (*p != '\n') ++p;
	}
	return hdr_from_dict(d);
}

bam_hdr_t *sam_hdr_read(htsFile *fp)
{
	if (fp->is_bin) {
		return bam_hdr_read(fp->fp.bgzf);
	} else if (fp->is_cram) {
		return cram_header_to_bam(fp->fp.cram->header);
	} else {
		kstring_t str;
		bam_hdr_t *h;
		int has_SQ = 0;
		str.l = str.m = 0; str.s = 0;
		while (hts_getline(fp, KS_SEP_LINE, &fp->line) >= 0) {
			if (fp->line.s[0] != '@') break;
			if (fp->line.l > 3 && strncmp(fp->line.s,"@SQ",3) == 0) has_SQ = 1;
			kputsn(fp->line.s, fp->line.l, &str);
			kputc('\n', &str);
		}
		if (! has_SQ && fp->fn_aux) {
			char line[2048];
			FILE *f = fopen(fp->fn_aux, "r");
			if (f == NULL) return NULL;
			while (fgets(line, sizeof line, f)) {
				const char *name = strtok(line, "\t");
				const char *length = strtok(NULL, "\t");
				ksprintf(&str, "@SQ\tSN:%s\tLN:%s\n", name, length);
			}
			fclose(f);
		}
		if (str.l == 0) kputsn("", 0, &str);
		h = sam_hdr_parse(str.l, str.s);
		h->l_text = str.l; h->text = str.s;
		return h;
	}
}

int sam_hdr_write(htsFile *fp, const bam_hdr_t *h)
{
	if (fp->is_bin) {
		bam_hdr_write(fp->fp.bgzf, h);
	} else if (fp->is_cram) {
		cram_fd *fd = fp->fp.cram;
		if (cram_set_header(fd, bam_header_to_cram((bam_hdr_t *)h)) < 0) return -1;
		if (cram_write_SAM_hdr(fd, fd->header) < 0) return -1;
	} else {
		char *p;
		hputs(h->text, fp->fp.hfile);
		p = strstr(h->text, "@SQ\t"); // FIXME: we need a loop to make sure "@SQ\t" does not match something unwanted!!!
		if (p == 0) {
			int i;
			for (i = 0; i < h->n_targets; ++i) {
				fp->line.l = 0;
				kputsn("@SQ\tSN:", 7, &fp->line); kputs(h->target_name[i], &fp->line);
				kputsn("\tLN:", 4, &fp->line); kputw(h->target_len[i], &fp->line); kputc('\n', &fp->line);
				hwrite(fp->fp.hfile, fp->line.s, fp->line.l);
			}
		}
		hflush(fp->fp.hfile);
	}
	return 0;
}

/**********************
 *** SAM record I/O ***
 **********************/

int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b)
{
#define _read_token(_p) (_p); for (; *(_p) && *(_p) != '\t'; ++(_p)); if (*(_p) != '\t') goto err_ret; *(_p)++ = 0
#define _read_token_aux(_p) (_p); for (; *(_p) && *(_p) != '\t'; ++(_p)); *(_p)++ = 0 // this is different in that it does not test *(_p)=='\t'
#define _get_mem(type_t, _x, _s, _l) ks_resize((_s), (_s)->l + (_l)); *(_x) = (type_t*)((_s)->s + (_s)->l); (_s)->l += (_l)
#define _parse_err(cond, msg) do { if ((cond) && hts_verbose >= 1) { fprintf(stderr, "[E::%s] " msg "\n", __func__); goto err_ret; } } while (0)
#define _parse_warn(cond, msg) if ((cond) && hts_verbose >= 2) fprintf(stderr, "[W::%s] " msg "\n", __func__)

	uint8_t *t;
	char *p = s->s, *q;
	int i;
	kstring_t str;
	bam1_core_t *c = &b->core;

	str.l = b->l_data = 0;
	str.s = (char*)b->data; str.m = b->m_data;
	memset(c, 0, 32);
	if (h->cigar_tab == 0) {
		h->cigar_tab = (uint8_t*)calloc(128, sizeof(uint8_t));
		for (i = 0; BAM_CIGAR_STR[i]; ++i)
			h->cigar_tab[(int)BAM_CIGAR_STR[i]] = i;
	}
	// qname
	q = _read_token(p);
	kputsn_(q, p - q, &str);
	c->l_qname = p - q;
	// flag
	c->flag = strtol(p, &p, 0);
	if (*p++ != '\t') goto err_ret; // malformated flag
	// chr
	q = _read_token(p);
	if (strcmp(q, "*")) {
		_parse_err(h->n_targets == 0, "missing SAM header");
		c->tid = bam_name2id(h, q);
		_parse_warn(c->tid < 0, "urecognized reference name; treated as unmapped");
	} else c->tid = -1;
	// pos
	c->pos = strtol(p, &p, 10) - 1;
	if (*p++ != '\t') goto err_ret;
	if (c->pos < 0 && c->tid >= 0) {
		_parse_warn(1, "mapped query cannot have zero coordinate; treated as unmapped");
		c->tid = -1;
	}
	if (c->tid < 0) c->flag |= BAM_FUNMAP;
	// mapq
	c->qual = strtol(p, &p, 10);
	if (*p++ != '\t') goto err_ret;
	// cigar
	if (*p != '*') {
		uint32_t *cigar;
		for (q = p, c->n_cigar = 0; *q && *q != '\t'; ++q)
			if (!isdigit(*q)) ++c->n_cigar;
		_get_mem(uint32_t, &cigar, &str, c->n_cigar<<2);
		for (i = 0, q = p; i < c->n_cigar; ++i, ++q) {
			int op;
			cigar[i] = strtol(q, &q, 10)<<BAM_CIGAR_SHIFT;
			op = (uint8_t)*q >= 128? -1 : h->cigar_tab[(int)*q];
			_parse_err(op < 0, "urecognized CIGAR operator");
			cigar[i] |= op;
		}
		p = q + 1;
		i = bam_cigar2rlen(c->n_cigar, cigar);
	} else {
		_parse_warn(!(c->flag&BAM_FUNMAP), "mapped query must have a CIGAR; treated as unmapped");
		c->flag |= BAM_FUNMAP;
		q = _read_token(p);
		i = 1;
	}
	c->bin = hts_reg2bin(c->pos, c->pos + i, 14, 5);
	// mate chr
	q = _read_token(p);
	if (strcmp(q, "=") == 0) c->mtid = c->tid;
	else if (strcmp(q, "*") == 0) c->mtid = -1;
	else c->mtid = bam_name2id(h, q);
	// mpos
	c->mpos = strtol(p, &p, 10) - 1;
	if (*p++ != '\t') goto err_ret;
	if (c->mpos < 0 && c->mtid >= 0) {
		_parse_warn(1, "mapped mate cannot have zero coordinate; treated as unmapped");
		c->mtid = -1;
	}
	// tlen
	c->isize = strtol(p, &p, 10);
	if (*p++ != '\t') goto err_ret;
	// seq
	q = _read_token(p);
	if (strcmp(q, "*")) {
		c->l_qseq = p - q - 1;
		i = bam_cigar2qlen(c->n_cigar, (uint32_t*)(str.s + c->l_qname));
		_parse_err(c->n_cigar && i != c->l_qseq, "CIGAR and query sequence are of different length");
		i = (c->l_qseq + 1) >> 1;
		_get_mem(uint8_t, &t, &str, i);
		memset(t, 0, i);
		for (i = 0; i < c->l_qseq; ++i)
			t[i>>1] |= seq_nt16_table[(int)q[i]] << ((~i&1)<<2);
	} else c->l_qseq = 0;
	// qual
	q = _read_token_aux(p);
	_get_mem(uint8_t, &t, &str, c->l_qseq);
	if (strcmp(q, "*")) {
		_parse_err(p - q - 1 != c->l_qseq, "SEQ and QUAL are of different length");
		for (i = 0; i < c->l_qseq; ++i) t[i] = q[i] - 33;
	} else memset(t, 0xff, c->l_qseq);
	// aux
	// Note that (like the bam1_core_t fields) this aux data in b->data is
	// stored in host endianness; so there is no byte swapping needed here.
	while (p < s->s + s->l) {
		uint8_t type;
		q = _read_token_aux(p); // FIXME: can be accelerated for long 'B' arrays
		_parse_err(p - q - 1 < 6, "incomplete aux field");
		kputsn_(q, 2, &str);
		q += 3; type = *q++; ++q; // q points to value
		if (type == 'A' || type == 'a' || type == 'c' || type == 'C') {
			kputc_('A', &str);
			kputc_(*q, &str);
		} else if (type == 'i' || type == 'I') {
			long x;
			x = strtol(q, &q, 10);
			if (x < 0) {
				if (x >= INT8_MIN) {
					kputc_('c', &str); kputc_(x, &str);
				} else if (x >= INT16_MIN) {
					int16_t y = x;
					kputc_('s', &str); kputsn_((char*)&y, 2, &str);
				} else {
					int32_t y = x;
					kputc_('i', &str); kputsn_(&y, 4, &str);
				}
			} else {
				if (x <= UINT8_MAX) {
					kputc_('C', &str); kputc_(x, &str);
				} else if (x <= UINT16_MAX) {
					uint16_t y = x;
					kputc_('S', &str); kputsn_(&y, 2, &str);
				} else {
					uint32_t y = x;
					kputc_('I', &str); kputsn_(&y, 4, &str);
				}
			}
		} else if (type == 'f') {
			float x;
			x = strtod(q, &q);
			kputc_('f', &str); kputsn_(&x, 4, &str);
		} else if (type == 'Z' || type == 'H') {
			kputc_(type, &str);kputsn_(q, p - q, &str); // note that this include the trailing NULL
		} else if (type == 'B') {
			int32_t n;
			char *r;
			_parse_err(p - q - 1 < 3, "incomplete B-typed aux field");
			type = *q++; // q points to the first ',' following the typing byte
			for (r = q, n = 0; *r; ++r)
				if (*r == ',') ++n;
			kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
			// FIXME: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
			if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0); kputc_(x, &str); }
			else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtoul(q + 1, &q, 0); kputc_(x, &str); }
			else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
			else if (type == 'S') while (q + 1 < p) { uint16_t x = strtoul(q + 1, &q, 0); kputsn_(&x, 2, &str); }
			else if (type == 'i') while (q + 1 < p) { int32_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
			else if (type == 'I') while (q + 1 < p) { uint32_t x = strtoul(q + 1, &q, 0); kputsn_(&x, 4, &str); }
			else if (type == 'f') while (q + 1 < p) { float    x = strtod(q + 1, &q);    kputsn_(&x, 4, &str); }
			else _parse_err(1, "unrecognized type");
		} else _parse_err(1, "unrecognized type");
	}
	b->data = (uint8_t*)str.s; b->l_data = str.l; b->m_data = str.m;
	return 0;

#undef _parse_warn
#undef _parse_err
#undef _get_mem
#undef _read_token_aux
#undef _read_token
err_ret:
	b->data = (uint8_t*)str.s; b->l_data = str.l; b->m_data = str.m;
	return -2;
}

int sam_read1(htsFile *fp, bam_hdr_t *h, bam1_t *b)
{
	if (fp->is_bin) {
		return bam_read1(fp->fp.bgzf, b);
	} else if (fp->is_cram) {
		return cram_get_bam_seq(fp->fp.cram, &b);
	} else {
		int ret;
err_recover:
		if (fp->line.l == 0) {
			ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
			if (ret < 0) return -1;
		}
		ret = sam_parse1(&fp->line, h, b);
		fp->line.l = 0;
		if (ret < 0) {
			if (hts_verbose >= 1)
				fprintf(stderr, "[W::%s] parse error at line %lld\n", __func__, (long long)fp->lineno);
			if (h->ignore_sam_err) goto err_recover;
		}
		return ret;
	}
}

int sam_format1(const bam_hdr_t *h, const bam1_t *b, kstring_t *str)
{
	int i;
	uint8_t *s;
	const bam1_core_t *c = &b->core;

	str->l = 0;
	kputsn(bam_get_qname(b), c->l_qname-1, str); kputc('\t', str); // query name
	kputw(c->flag, str); kputc('\t', str); // flag
	if (c->tid >= 0) { // chr
		kputs(h->target_name[c->tid] , str);
		kputc('\t', str);
	} else kputsn("*\t", 2, str);
	kputw(c->pos + 1, str); kputc('\t', str); // pos
	kputw(c->qual, str); kputc('\t', str); // qual
	if (c->n_cigar) { // cigar
		uint32_t *cigar = bam_get_cigar(b);
		for (i = 0; i < c->n_cigar; ++i) {
			kputw(bam_cigar_oplen(cigar[i]), str);
			kputc(bam_cigar_opchr(cigar[i]), str);
		}
	} else kputc('*', str);
	kputc('\t', str);
	if (c->mtid < 0) kputsn("*\t", 2, str); // mate chr
	else if (c->mtid == c->tid) kputsn("=\t", 2, str);
	else {
		kputs(h->target_name[c->mtid], str);
		kputc('\t', str);
	}
	kputw(c->mpos + 1, str); kputc('\t', str); // mate pos
	kputw(c->isize, str); kputc('\t', str); // template len
	if (c->l_qseq) { // seq and qual
		uint8_t *s = bam_get_seq(b);
		for (i = 0; i < c->l_qseq; ++i) kputc("=ACMGRSVTWYHKDBN"[bam_seqi(s, i)], str);
		kputc('\t', str);
		s = bam_get_qual(b);
		if (s[0] == 0xff) kputc('*', str);
		else for (i = 0; i < c->l_qseq; ++i) kputc(s[i] + 33, str);
	} else kputsn("*\t*", 3, str);
	s = bam_get_aux(b); // aux
	while (s < b->data + b->l_data) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s++;
		kputc('\t', str); kputsn((char*)key, 2, str); kputc(':', str);
		if (type == 'A') { kputsn("A:", 2, str); kputc(*s, str); ++s; }
		else if (type == 'C') { kputsn("i:", 2, str); kputw(*s, str); ++s; }
		else if (type == 'c') { kputsn("i:", 2, str); kputw(*(int8_t*)s, str); ++s; }
		else if (type == 'S') { kputsn("i:", 2, str); kputw(*(uint16_t*)s, str); s += 2; }
		else if (type == 's') { kputsn("i:", 2, str); kputw(*(int16_t*)s, str); s += 2; }
		else if (type == 'I') { kputsn("i:", 2, str); kputuw(*(uint32_t*)s, str); s += 4; }
		else if (type == 'i') { kputsn("i:", 2, str); kputw(*(int32_t*)s, str); s += 4; }
		else if (type == 'f') { ksprintf(str, "f:%g", *(float*)s); s += 4; }
		else if (type == 'd') { ksprintf(str, "d:%g", *(double*)s); s += 8; }
		else if (type == 'Z' || type == 'H') { kputc(type, str); kputc(':', str); while (*s) kputc(*s++, str); ++s; }
		else if (type == 'B') {
			uint8_t sub_type = *(s++);
			int32_t n;
			memcpy(&n, s, 4);
			s += 4; // no point to the start of the array
			kputsn("B:", 2, str); kputc(sub_type, str); // write the typing
			for (i = 0; i < n; ++i) { // FIXME: for better performance, put the loop after "if"
				kputc(',', str);
				if ('c' == sub_type)      { kputw(*(int8_t*)s, str); ++s; }
				else if ('C' == sub_type) { kputw(*(uint8_t*)s, str); ++s; }
				else if ('s' == sub_type) { kputw(*(int16_t*)s, str); s += 2; }
				else if ('S' == sub_type) { kputw(*(uint16_t*)s, str); s += 2; }
				else if ('i' == sub_type) { kputw(*(int32_t*)s, str); s += 4; }
				else if ('I' == sub_type) { kputuw(*(uint32_t*)s, str); s += 4; }
				else if ('f' == sub_type) { ksprintf(str, "%g", *(float*)s); s += 4; }
			}
		}
	}
	return str->l;
}

int sam_write1(htsFile *fp, const bam_hdr_t *h, const bam1_t *b)
{
	if (fp->is_bin) {
		return bam_write1(fp->fp.bgzf, b);
	} else if (fp->is_cram) {
		return cram_put_bam_seq(fp->fp.cram, (bam1_t *)b);
	} else {
		sam_format1(h, b, &fp->line);
		hwrite(fp->fp.hfile, fp->line.s, fp->line.l);
		hputc('\n', fp->fp.hfile);
		return fp->line.l + 1;
	}
}

/************************
 *** Auxiliary fields ***
 ************************/

int bam_aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f') return 4;
	else return 0;
}

void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, uint8_t *data)
{
	int ori_len = b->l_data;
	b->l_data += 3 + len;
	if (b->m_data < b->l_data) {
		b->m_data = b->l_data;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	b->data[ori_len] = tag[0]; b->data[ori_len + 1] = tag[1];
	b->data[ori_len + 2] = type;
	memcpy(b->data + ori_len + 3, data, len);
}

#define __skip_tag(s) do { \
		int type = toupper(*(s)); \
		++(s); \
		if (type == 'Z' || type == 'H') { while (*(s)) ++(s); ++(s); } \
		else if (type == 'B') (s) += 5 + bam_aux_type2size(*(s)) * (*(int32_t*)((s)+1)); \
		else (s) += bam_aux_type2size(type); \
	} while(0)

uint8_t *bam_aux_get(const bam1_t *b, const char tag[2])
{
	uint8_t *s;
	int y = tag[0]<<8 | tag[1];
	s = bam_get_aux(b);
	while (s < b->data + b->l_data) {
		int x = (int)s[0]<<8 | s[1];
		s += 2;
		if (x == y) return s;
		__skip_tag(s);
	}
	return 0;
}
// s MUST BE returned by bam_aux_get()
int bam_aux_del(bam1_t *b, uint8_t *s)
{
	uint8_t *p, *aux;
	int l_aux = bam_get_l_aux(b);
	aux = bam_get_aux(b);
	p = s - 2;
	__skip_tag(s);
	memmove(p, s, l_aux - (s - aux));
	b->l_data -= s - p;
	return 0;
}

int32_t bam_aux2i(const uint8_t *s)
{
	int type;
	type = *s++;
	if (type == 'c') return (int32_t)*(int8_t*)s;
	else if (type == 'C') return (int32_t)*(uint8_t*)s;
	else if (type == 's') return (int32_t)*(int16_t*)s;
	else if (type == 'S') return (int32_t)*(uint16_t*)s;
	else if (type == 'i' || type == 'I') return *(int32_t*)s;
	else return 0;
}

double bam_aux2f(const uint8_t *s)
{
	int type;
	type = *s++;
	if (type == 'd') return *(double*)s;
	else if (type == 'f') return *(float*)s;
	else return 0.0;
}

char bam_aux2A(const uint8_t *s)
{
	int type;
	type = *s++;
	if (type == 'A') return *(char*)s;
	else return 0;
}

char *bam_aux2Z(const uint8_t *s)
{
	int type;
	type = *s++;
	if (type == 'Z' || type == 'H') return (char*)s;
	else return 0;
}

/**************************
 *** Pileup and Mpileup ***
 **************************/

#if !defined(BAM_NO_PILEUP)

#include <assert.h>

#define BAM_DEF_MASK (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

/*******************
 *** Memory pool ***
 *******************/

typedef struct {
	int k, x, y, end;
} cstate_t;

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

typedef struct __linkbuf_t {
	bam1_t b;
	int32_t beg, end;
	cstate_t s;
	struct __linkbuf_t *next;
} lbnode_t;

typedef struct {
	int cnt, n, max;
	lbnode_t **buf;
} mempool_t;

static mempool_t *mp_init()
{
	mempool_t *mp;
	mp = (mempool_t*)calloc(1, sizeof(mempool_t));
	return mp;
}
static void mp_destroy(mempool_t *mp)
{
	int k;
	for (k = 0; k < mp->n; ++k) {
		free(mp->buf[k]->b.data);
		free(mp->buf[k]);
	}
	free(mp->buf);
	free(mp);
}
static inline lbnode_t *mp_alloc(mempool_t *mp)
{
	++mp->cnt;
	if (mp->n == 0) return (lbnode_t*)calloc(1, sizeof(lbnode_t));
	else return mp->buf[--mp->n];
}
static inline void mp_free(mempool_t *mp, lbnode_t *p)
{
	--mp->cnt; p->next = 0; // clear lbnode_t::next here
	if (mp->n == mp->max) {
		mp->max = mp->max? mp->max<<1 : 256;
		mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
	}
	mp->buf[mp->n++] = p;
}

/**********************
 *** CIGAR resolver ***
 **********************/

/* s->k: the index of the CIGAR operator that has just been processed.
   s->x: the reference coordinate of the start of s->k
   s->y: the query coordiante of the start of s->k
 */
static inline int resolve_cigar2(bam_pileup1_t *p, int32_t pos, cstate_t *s)
{
#define _cop(c) ((c)&BAM_CIGAR_MASK)
#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)

	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t *cigar = bam_get_cigar(b);
	int k;
	// determine the current CIGAR operation
//	fprintf(stderr, "%s\tpos=%d\tend=%d\t(%d,%d,%d)\n", bam1_qname(b), pos, s->end, s->k, s->x, s->y);
	if (s->k == -1) { // never processed
		if (c->n_cigar == 1) { // just one operation, save a loop
		  if (_cop(cigar[0]) == BAM_CMATCH || _cop(cigar[0]) == BAM_CEQUAL || _cop(cigar[0]) == BAM_CDIFF) s->k = 0, s->x = c->pos, s->y = 0;
		} else { // find the first match or deletion
			for (k = 0, s->x = c->pos, s->y = 0; k < c->n_cigar; ++k) {
				int op = _cop(cigar[k]);
				int l = _cln(cigar[k]);
				if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CEQUAL || op == BAM_CDIFF) break;
				else if (op == BAM_CREF_SKIP) s->x += l;
				else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
			}
			assert(k < c->n_cigar);
			s->k = k;
		}
	} else { // the read has been processed before
		int op, l = _cln(cigar[s->k]);
		if (pos - s->x >= l) { // jump to the next operation
			assert(s->k < c->n_cigar); // otherwise a bug: this function should not be called in this case
			op = _cop(cigar[s->k+1]);
			if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) { // jump to the next without a loop
			  if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
				s->x += l;
				++s->k;
			} else { // find the next M/D/N/=/X
			  if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
				s->x += l;
				for (k = s->k + 1; k < c->n_cigar; ++k) {
					op = _cop(cigar[k]), l = _cln(cigar[k]);
					if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) break;
					else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
				}
				s->k = k;
			}
			assert(s->k < c->n_cigar); // otherwise a bug
		} // else, do nothing
	}
	{ // collect pileup information
		int op, l;
		op = _cop(cigar[s->k]); l = _cln(cigar[s->k]);
		p->is_del = p->indel = p->is_refskip = 0;
		if (s->x + l - 1 == pos && s->k + 1 < c->n_cigar) { // peek the next operation
			int op2 = _cop(cigar[s->k+1]);
			int l2 = _cln(cigar[s->k+1]);
			if (op2 == BAM_CDEL) p->indel = -(int)l2;
			else if (op2 == BAM_CINS) p->indel = l2;
			else if (op2 == BAM_CPAD && s->k + 2 < c->n_cigar) { // no working for adjacent padding
				int l3 = 0;
				for (k = s->k + 2; k < c->n_cigar; ++k) {
					op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
					if (op2 == BAM_CINS) l3 += l2;
					else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF) break;
				}
				if (l3 > 0) p->indel = l3;
			}
		}
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			p->qpos = s->y + (pos - s->x);
		} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
			p->is_del = 1; p->qpos = s->y; // FIXME: distinguish D and N!!!!!
			p->is_refskip = (op == BAM_CREF_SKIP);
		} // cannot be other operations; otherwise a bug
		p->is_head = (pos == c->pos); p->is_tail = (pos == s->end);
	}
	return 1;
}

/***********************
 *** Pileup iterator ***
 ***********************/

struct __bam_plp_t {
	mempool_t *mp;
	lbnode_t *head, *tail, *dummy;
	int32_t tid, pos, max_tid, max_pos;
	int is_eof, flag_mask, max_plp, error, maxcnt;
	uint64_t id;
	bam_pileup1_t *plp;
	// for the "auto" interface only
	bam1_t *b;
	bam_plp_auto_f func;
	void *data;
};

bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data)
{
	bam_plp_t iter;
	iter = (bam_plp_t)calloc(1, sizeof(struct __bam_plp_t));
	iter->mp = mp_init();
	iter->head = iter->tail = mp_alloc(iter->mp);
	iter->dummy = mp_alloc(iter->mp);
	iter->max_tid = iter->max_pos = -1;
	iter->flag_mask = BAM_DEF_MASK;
	iter->maxcnt = 8000;
	if (func) {
		iter->func = func;
		iter->data = data;
		iter->b = bam_init1();
	}
	return iter;
}

void bam_plp_destroy(bam_plp_t iter)
{
	mp_free(iter->mp, iter->dummy);
	mp_free(iter->mp, iter->head);
	if (iter->mp->cnt != 0)
		fprintf(stderr, "[bam_plp_destroy] memory leak: %d. Continue anyway.\n", iter->mp->cnt);
	mp_destroy(iter->mp);
	if (iter->b) bam_destroy1(iter->b);
	free(iter->plp);
	free(iter);
}

// Prepares next pileup position in bam records collected by bam_plp_auto -> user func -> bam_plp_push. Returns
// pointer to the piled records if next position is ready or NULL if there is not enough records in the
// buffer yet (the current position is still the maximum position across all buffered reads).
const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
	if (iter->error) { *_n_plp = -1; return 0; }
	*_n_plp = 0;
	if (iter->is_eof && iter->head->next == 0) return 0;
	while (iter->is_eof || iter->max_tid > iter->tid || (iter->max_tid == iter->tid && iter->max_pos > iter->pos)) {
		int n_plp = 0;
		lbnode_t *p, *q;
		// write iter->plp at iter->pos
		iter->dummy->next = iter->head;
		for (p = iter->head, q = iter->dummy; p->next; q = p, p = p->next) {
			if (p->b.core.tid < iter->tid || (p->b.core.tid == iter->tid && p->end <= iter->pos)) { // then remove
				q->next = p->next; mp_free(iter->mp, p); p = q;
			} else if (p->b.core.tid == iter->tid && p->beg <= iter->pos) { // here: p->end > pos; then add to pileup
				if (n_plp == iter->max_plp) { // then double the capacity
					iter->max_plp = iter->max_plp? iter->max_plp<<1 : 256;
					iter->plp = (bam_pileup1_t*)realloc(iter->plp, sizeof(bam_pileup1_t) * iter->max_plp);
				}
				iter->plp[n_plp].b = &p->b;
				if (resolve_cigar2(iter->plp + n_plp, iter->pos, &p->s)) ++n_plp; // actually always true...
			}
		}
		iter->head = iter->dummy->next; // dummy->next may be changed
		*_n_plp = n_plp; *_tid = iter->tid; *_pos = iter->pos;
		// update iter->tid and iter->pos
		if (iter->head->next) {
			if (iter->tid > iter->head->b.core.tid) {
				fprintf(stderr, "[%s] unsorted input. Pileup aborts.\n", __func__);
				iter->error = 1;
				*_n_plp = -1;
				return 0;
			}
		}
		if (iter->tid < iter->head->b.core.tid) { // come to a new reference sequence
			iter->tid = iter->head->b.core.tid; iter->pos = iter->head->beg; // jump to the next reference
		} else if (iter->pos < iter->head->beg) { // here: tid == head->b.core.tid
			iter->pos = iter->head->beg; // jump to the next position
		} else ++iter->pos; // scan contiguously
		// return
		if (n_plp) return iter->plp;
		if (iter->is_eof && iter->head->next == 0) break;
	}
	return 0;
}

int bam_plp_push(bam_plp_t iter, const bam1_t *b)
{
	if (iter->error) return -1;
	if (b) {
		if (b->core.tid < 0) return 0;
		if (b->core.flag & iter->flag_mask) return 0;
		if (iter->tid == b->core.tid && iter->pos == b->core.pos && iter->mp->cnt > iter->maxcnt) return 0;
		bam_copy1(&iter->tail->b, b);
#ifndef BAM_NO_ID
		iter->tail->b.id = iter->id++;
#endif
		iter->tail->beg = b->core.pos;
		iter->tail->end = b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
		iter->tail->s = g_cstate_null; iter->tail->s.end = iter->tail->end - 1; // initialize cstate_t
		if (b->core.tid < iter->max_tid) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (chromosomes out of order)\n");
			iter->error = 1;
			return -1;
		}
		if ((b->core.tid == iter->max_tid) && (iter->tail->beg < iter->max_pos)) {
			fprintf(stderr, "[bam_pileup_core] the input is not sorted (reads out of order)\n");
			iter->error = 1;
			return -1;
		}
		iter->max_tid = b->core.tid; iter->max_pos = iter->tail->beg;
		if (iter->tail->end > iter->pos || iter->tail->b.core.tid > iter->tid) {
			iter->tail->next = mp_alloc(iter->mp);
			iter->tail = iter->tail->next;
		}
	} else iter->is_eof = 1;
	return 0;
}

const bam_pileup1_t *bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
	const bam_pileup1_t *plp;
	if (iter->func == 0 || iter->error) { *_n_plp = -1; return 0; }
	if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
	else { // no pileup line can be obtained; read alignments
		*_n_plp = 0;
		if (iter->is_eof) return 0;
		while (iter->func(iter->data, iter->b) >= 0) {
			if (bam_plp_push(iter, iter->b) < 0) {
				*_n_plp = -1;
				return 0;
			}
			if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
			// otherwise no pileup line can be returned; read the next alignment.
		}
		bam_plp_push(iter, 0);
		if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
		return 0;
	}
}

void bam_plp_reset(bam_plp_t iter)
{
	lbnode_t *p, *q;
	iter->max_tid = iter->max_pos = -1;
	iter->tid = iter->pos = 0;
	iter->is_eof = 0;
	for (p = iter->head; p->next;) {
		q = p->next;
		mp_free(iter->mp, p);
		p = q;
	}
	iter->head = iter->tail;
}

void bam_plp_set_mask(bam_plp_t iter, int mask)
{
	iter->flag_mask = mask < 0? BAM_DEF_MASK : (BAM_FUNMAP | mask);
}

void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt)
{
	iter->maxcnt = maxcnt;
}

/************************
 *** Mpileup iterator ***
 ************************/

struct __bam_mplp_t {
	int n;
	uint64_t min, *pos;
	bam_plp_t *iter;
	int *n_plp;
	const bam_pileup1_t **plp;
};

bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data)
{
	int i;
	bam_mplp_t iter;
	iter = (bam_mplp_t)calloc(1, sizeof(struct __bam_mplp_t));
	iter->pos = (uint64_t*)calloc(n, 8);
	iter->n_plp = (int*)calloc(n, sizeof(int));
	iter->plp = (const bam_pileup1_t**)calloc(n, sizeof(void*));
	iter->iter = (bam_plp_t*)calloc(n, sizeof(void*));
	iter->n = n;
	iter->min = (uint64_t)-1;
	for (i = 0; i < n; ++i) {
		iter->iter[i] = bam_plp_init(func, data[i]);
		iter->pos[i] = iter->min;
	}
	return iter;
}

void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt)
{
	int i;
	for (i = 0; i < iter->n; ++i)
		iter->iter[i]->maxcnt = maxcnt;
}

void bam_mplp_destroy(bam_mplp_t iter)
{
	int i;
	for (i = 0; i < iter->n; ++i) bam_plp_destroy(iter->iter[i]);
	free(iter->iter); free(iter->pos); free(iter->n_plp); free(iter->plp);
	free(iter);
}

int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp)
{
	int i, ret = 0;
	uint64_t new_min = (uint64_t)-1;
	for (i = 0; i < iter->n; ++i) {
		if (iter->pos[i] == iter->min) {
			int tid, pos;
			iter->plp[i] = bam_plp_auto(iter->iter[i], &tid, &pos, &iter->n_plp[i]);
			iter->pos[i] = iter->plp[i] ? (uint64_t)tid<<32 | pos : 0;
		}
		if (iter->plp[i] && iter->pos[i] < new_min) new_min = iter->pos[i];
	}
	iter->min = new_min;
	if (new_min == (uint64_t)-1) return 0;
	*_tid = new_min>>32; *_pos = (uint32_t)new_min;
	for (i = 0; i < iter->n; ++i) {
		if (iter->pos[i] == iter->min) { // FIXME: valgrind reports "uninitialised value(s) at this line"
			n_plp[i] = iter->n_plp[i], plp[i] = iter->plp[i];
			++ret;
		} else n_plp[i] = 0, plp[i] = 0;
	}
	return ret;
}

#endif // ~!defined(BAM_NO_PILEUP)
