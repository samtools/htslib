HTS_DIR = ../
include $(HTS_DIR)/htslib_static.mk

CC = gcc
CFLAGS = -Wall -O2

#to statically link to libhts
LDFLAGS = $(HTS_DIR)/libhts.a -L$(HTS_DIR) $(HTSLIB_static_LDFLAGS) $(HTSLIB_static_LIBS)

#to dynamically link to libhts
#LDFLAGS = -L $(HTS_DIR) -lhts -Wl,-rpath,<installed lib path>

PRGS = flags split split2 cram read_fast read_header read_ref read_bam \
	read_aux dump_aux add_header rem_header update_header mod_bam mod_aux \
	mod_aux_ba write_fast idx_on_write read_reg read_multireg pileup \
	mpileup modstate pileup_mod flags_field split_t1 split_t2 \
	read_fast_i qtask_ordered qtask_unordered index_fasta

all: $(PRGS)

flags: flags_demo.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) flags_demo.c -o $@  $(LDFLAGS)

split: split.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) split.c -o $@  $(LDFLAGS)

split2: split2.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) split2.c -o $@  $(LDFLAGS)

cram: cram.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) cram.c -o $@  $(LDFLAGS)

read_fast: read_fast.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_fast.c -o $@  $(LDFLAGS)

read_header: read_header.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_header.c -o $@  $(LDFLAGS)

read_ref: read_refname.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_refname.c -o $@  $(LDFLAGS)

read_bam: read_bam.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_bam.c -o $@  $(LDFLAGS)

read_aux: read_aux.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_aux.c -o $@  $(LDFLAGS)

dump_aux: dump_aux.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) dump_aux.c -o $@  $(LDFLAGS)

add_header: add_header.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) add_header.c -o $@  $(LDFLAGS)

rem_header: rem_header.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) rem_header.c -o $@  $(LDFLAGS)

update_header: update_header.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) update_header.c -o $@  $(LDFLAGS)

mod_bam: mod_bam.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) mod_bam.c -o $@  $(LDFLAGS)

mod_aux: mod_aux.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) mod_aux.c -o $@  $(LDFLAGS)

mod_aux_ba: mod_aux_ba.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) mod_aux_ba.c -o $@  $(LDFLAGS)

write_fast: write_fast.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) write_fast.c -o $@  $(LDFLAGS)

idx_on_write: index_write.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) index_write.c -o $@  $(LDFLAGS)

read_reg: index_reg_read.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) index_reg_read.c -o $@  $(LDFLAGS)

read_multireg: index_multireg_read.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) index_multireg_read.c -o $@  $(LDFLAGS)

read_fast_i: read_fast_index.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_fast_index.c -o $@  $(LDFLAGS)

pileup: pileup.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) pileup.c -o $@  $(LDFLAGS)

mpileup: mpileup.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) mpileup.c -o $@  $(LDFLAGS)

modstate: modstate.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) modstate.c -o $@  $(LDFLAGS)

pileup_mod: pileup_mod.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) pileup_mod.c -o $@  $(LDFLAGS)

flags_field: flags_htsopt_field.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) flags_htsopt_field.c -o $@  $(LDFLAGS)

split_t1: split_thread1.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) split_thread1.c -o $@  $(LDFLAGS)

split_t2: split_thread2.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) split_thread2.c -o $@  $(LDFLAGS)

index_fasta: index_fasta.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) index_fasta.c -o $@  $(LDFLAGS)

qtask_ordered: qtask_ordered.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) qtask_ordered.c -o $@  $(LDFLAGS)

qtask_unordered: qtask_unordered.c
	$(CC) $(CFLAGS) -I $(HTS_DIR) qtask_unordered.c -o $@  $(LDFLAGS)

clean:
	find . -name "*.o" | xargs rm -rf
	find . -name "*.dSYM" | xargs rm -rf
	-rm -f $(PRGS)
