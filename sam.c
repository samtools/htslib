#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <zlib.h>
#include "endian.h"
#include "bgzf.h"
#include "kstring.h"
#include "sam.h"

#include "kseq.h"
KSTREAM_DECLARE(gzFile, gzread)

#include "khash.h"
KHASH_DECLARE(s2i, kh_cstr_t, int)

typedef khash_t(s2i) sdict_t;

/**************
 * Header I/O *
 **************/

sam_hdr_t *sam_hdr_init()
{
	return (sam_hdr_t*)calloc(1, sizeof(sam_hdr_t));
}

void sam_hdr_destroy(sam_hdr_t *h)
{
	int32_t i;
	if (h == 0) return;
	if (h->target_name) {
		for (i = 0; i < h->n_targets; ++i)
			free(h->target_name[i]);
		free(h->target_name);
		free(h->target_len);
	}
	free(h->text);
	if (h->dict) kh_destroy(s2i, (sdict_t*)h->dict);
	free(h);
}

sam_hdr_t *sam_hdr_read(htsFile *fp)
{
	sam_hdr_t *h;
	if (fp->is_bin) {
		char buf[4];
		int magic_len, has_EOF;
		int32_t i = 1, name_len;
		// check EOF
		has_EOF = bgzf_check_EOF((BGZF*)fp->fp);
		if (has_EOF < 0) {
			if (errno != ESPIPE && hts_verbose >= 2) perror("[W::sam_hdr_read] bgzf_check_EOF");
		} else if (has_EOF == 0 && hts_verbose >= 2)
			fprintf(stderr, "[W::%s] EOF marker is absent. The input is probably truncated.\n", __func__);
		// read "BAM1"
		magic_len = bgzf_read((BGZF*)fp->fp, buf, 4);
		if (magic_len != 4 || strncmp(buf, "BAM\1", 4)) {
			if (hts_verbose >= 1) fprintf(stderr, "[E::%s] invalid BAM binary header\n", __func__);
			return 0;
		}
		h = sam_hdr_init();
		// read plain text and the number of reference sequences
		bgzf_read((BGZF*)fp->fp, &h->l_text, 4);
		if (fp->is_be) ed_swap_4p(&h->l_text);
		h->text = (char*)calloc(h->l_text + 1, 1);
		bgzf_read((BGZF*)fp->fp, h->text, h->l_text);
		bgzf_read((BGZF*)fp->fp, &h->n_targets, 4);
		if (fp->is_be) ed_swap_4p(&h->n_targets);
		// read reference sequence names and lengths
		h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
		h->target_len = (uint32_t*)calloc(h->n_targets, 4);
		for (i = 0; i != h->n_targets; ++i) {
			bgzf_read((BGZF*)fp->fp, &name_len, 4);
			if (fp->is_be) ed_swap_4p(&name_len);
			h->target_name[i] = (char*)calloc(name_len, 1);
			bgzf_read((BGZF*)fp->fp, h->target_name[i], name_len);
			bgzf_read((BGZF*)fp->fp, &h->target_len[i], 4);
			if (fp->is_be) ed_swap_4p(&h->target_len[i]);
		}
		return h;
	}
	return 0;
}

int sam_hdr_write(htsFile *fp, const sam_hdr_t *h)
{
	if (fp->is_bin) {
		char buf[4];
		int32_t i, name_len, x;
		// write "BAM1"
		strncpy(buf, "BAM\1", 4);
		bgzf_write((BGZF*)fp->fp, buf, 4);
		// write plain text and the number of reference sequences
		if (fp->is_be) {
			x = ed_swap_4(h->l_text);
			bgzf_write((BGZF*)fp->fp, &x, 4);
			if (h->l_text) bgzf_write((BGZF*)fp->fp, h->text, h->l_text);
			x = ed_swap_4(h->n_targets);
			bgzf_write((BGZF*)fp->fp, &x, 4);
		} else {
			bgzf_write((BGZF*)fp->fp, &h->l_text, 4);
			if (h->l_text) bgzf_write((BGZF*)fp->fp, h->text, h->l_text);
			bgzf_write((BGZF*)fp->fp, &h->n_targets, 4);
		}
		// write sequence names and lengths
		for (i = 0; i != h->n_targets; ++i) {
			char *p = h->target_name[i];
			name_len = strlen(p) + 1;
			if (fp->is_be) {
				x = ed_swap_4(name_len);
				bgzf_write((BGZF*)fp->fp, &x, 4);
			} else bgzf_write((BGZF*)fp->fp, &name_len, 4);
			bgzf_write((BGZF*)fp->fp, p, name_len);
			if (fp->is_be) {
				x = ed_swap_4(h->target_len[i]);
				bgzf_write((BGZF*)fp->fp, &x, 4);
			} else bgzf_write((BGZF*)fp->fp, &h->target_len[i], 4);
		}
		bgzf_flush((BGZF*)fp->fp);
	}
	return 0;
}

/******************
 * SAM record I/O *
 ******************/

sam1_t *sam_init1()
{
	return (sam1_t*)calloc(1, sizeof(sam1_t));
}

static inline int sam_aux_type2size(int x)
{
	if (x == 'C' || x == 'c' || x == 'A') return 1;
	else if (x == 'S' || x == 's') return 2;
	else if (x == 'I' || x == 'i' || x == 'f') return 4;
	else return 0;
}

static void swap_data(const sam1_core_t *c, int l_data, uint8_t *data)
{
	uint8_t *s;
	uint32_t i, *cigar = (uint32_t*)(data + c->l_qname);
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
			int32_t n, Bsize = sam_aux_type2size(*s);
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

int sam_read1(htsFile *fp, const sam_hdr_t *h, sam1_t *b)
{
	if (fp->is_bin) {
		sam1_core_t *c = &b->core;
		int32_t block_len, ret, i;
		uint32_t x[8];

		if ((ret = bgzf_read((BGZF*)fp->fp, &block_len, 4)) != 4) {
			if (ret == 0) return -1; // normal end-of-file
			else return -2; // truncated
		}
		if (bgzf_read((BGZF*)fp->fp, x, 32) != 32) return -3;
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
		if (bgzf_read((BGZF*)fp->fp, b->data, b->l_data) != b->l_data) return -4;
		b->l_aux = b->l_data - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq+1)/2;
		if (fp->is_be) swap_data(c, b->l_data, b->data);
		return 4 + block_len;
	} else return 0;
}

int sam_format1(const sam_hdr_t *h, const sam1_t *b, kstring_t *str)
{
	int i;
	uint8_t *s;
	const sam1_core_t *c = &b->core;

	str->l = 0;
	kputsn(sam_get_qname(b), c->l_qname-1, str); kputc('\t', str); // query name
	kputw(c->flag, str); kputc('\t', str); // flag
	if (c->tid >= 0) { // chr
		kputs(h->target_name[c->tid] , str);
		kputc('\t', str);
	} else kputsn("*\t", 2, str);
	kputw(c->pos + 1, str); kputc('\t', str); // pos
	kputw(c->qual, str); kputc('\t', str); // qual
	if (c->n_cigar) { // cigar
		uint32_t *cigar = sam_get_cigar(b);
		for (i = 0; i < c->n_cigar; ++i) {
			kputw(sam_cigar_oplen(cigar[i]), str);
			kputc(sam_cigar_opchr(cigar[i]), str);
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
		uint8_t *s = sam_get_seq(b);
		for (i = 0; i < c->l_qseq; ++i) kputc("=ACMGRSVTWYHKDBN"[sam_seqi(s, i)], str);
		kputc('\t', str);
		s = sam_get_qual(b);
		if (s[0] == 0xff) kputc('*', str);
		else for (i = 0; i < c->l_qseq; ++i) kputc(s[i] + 33, str);
	} else kputsn("*\t*", 3, str);
	s = sam_get_aux(b); // aux
	while (s < b->data + b->l_data) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s; ++s;
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
			kputc(type, str); kputc(':', str); kputc(sub_type, str); // write the typing
			for (i = 0; i < n; ++i) {
				kputc(',', str);
				if ('c' == sub_type || 'c' == sub_type) { kputw(*(int8_t*)s, str); ++s; }
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

int sam_write1(htsFile *fp, const sam_hdr_t *h, const sam1_t *b)
{
	if (fp->is_bin) {
		const sam1_core_t *c = &b->core;
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
		bgzf_flush_try((BGZF*)fp->fp, 4 + block_len);
		if (fp->is_be) {
			for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
			y = block_len;
			bgzf_write((BGZF*)fp->fp, ed_swap_4p(&y), 4);
			swap_data(c, b->l_data, b->data);
		} else bgzf_write((BGZF*)fp->fp, &block_len, 4);
		bgzf_write((BGZF*)fp->fp, x, 32);
		bgzf_write((BGZF*)fp->fp, b->data, b->l_data);
		if (fp->is_be) swap_data(c, b->l_data, b->data);
		return 4 + block_len;
	} else {
		sam_format1(h, b, &fp->line);
		fwrite(fp->line.s, 1, fp->line.l, (FILE*)fp->fp);
		fputc('\n', (FILE*)fp->fp);
		return fp->line.l + 1;
	}
}
