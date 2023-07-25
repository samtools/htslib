HTS_DIR = ../
include $(HTS_DIR)/htslib_static.mk

CC = gcc
CFLAGS = -Wall -g -O0

#to statically link to libhts
LDFLAGS = $(HTS_DIR)/libhts.a -L$(HTS_DIR) $(HTSLIB_static_LDFLAGS) $(HTSLIB_static_LIBS)

#to dynamically link to libhts
#LDFLAGS = -L $(HTS_DIR) -lhts -Wl,-rpath,<installed lib path>

PRGS = flags split split2 cram read_fast read_header read_ref read_bam \
	read_aux dump_aux add_header rem_header update_header mod_bam mod_aux \
	mod_aux_ba write_fast idx_on_write read_reg read_multireg pileup \
	mpileup modstate pileup_mod flags_field split_t1 split_t2

all: $(PRGS)

flags:
	$(CC) $(CFLAGS) -I $(HTS_DIR) flags_demo.c -o $@  $(LDFLAGS)

split:
	$(CC) $(CFLAGS) -I $(HTS_DIR) split.c -o $@  $(LDFLAGS)

split2:
	$(CC) $(CFLAGS) -I $(HTS_DIR) split2.c -o $@  $(LDFLAGS)

cram:
	$(CC) $(CFLAGS) -I $(HTS_DIR) cram.c -o $@  $(LDFLAGS)

read_fast:
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_fast.c -o $@  $(LDFLAGS)

read_header:
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_header.c -o $@  $(LDFLAGS)

read_ref:
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_refname.c -o $@  $(LDFLAGS)

read_bam:
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_bam.c -o $@  $(LDFLAGS)

read_aux:
	$(CC) $(CFLAGS) -I $(HTS_DIR) read_aux.c -o $@  $(LDFLAGS)

dump_aux:
	$(CC) $(CFLAGS) -I $(HTS_DIR) dump_aux.c -o $@  $(LDFLAGS)

add_header:
	$(CC) $(CFLAGS) -I $(HTS_DIR) add_header.c -o $@  $(LDFLAGS)

rem_header:
	$(CC) $(CFLAGS) -I $(HTS_DIR) rem_header.c -o $@  $(LDFLAGS)

update_header:
	$(CC) $(CFLAGS) -I $(HTS_DIR) update_header.c -o $@  $(LDFLAGS)

mod_bam:
	$(CC) $(CFLAGS) -I $(HTS_DIR) mod_bam.c -o $@  $(LDFLAGS)

mod_aux:
	$(CC) $(CFLAGS) -I $(HTS_DIR) mod_aux.c -o $@  $(LDFLAGS)

mod_aux_ba:
	$(CC) $(CFLAGS) -I $(HTS_DIR) mod_aux_ba.c -o $@  $(LDFLAGS)

write_fast:
	$(CC) $(CFLAGS) -I $(HTS_DIR) write_fast.c -o $@  $(LDFLAGS)

idx_on_write:
	$(CC) $(CFLAGS) -I $(HTS_DIR) index_write.c -o $@  $(LDFLAGS)

read_reg:
	$(CC) $(CFLAGS) -I $(HTS_DIR) index_reg_read.c -o $@  $(LDFLAGS)

read_multireg:
	$(CC) $(CFLAGS) -I $(HTS_DIR) index_multireg_read.c -o $@  $(LDFLAGS)

pileup:
	$(CC) $(CFLAGS) -I $(HTS_DIR) pileup.c -o $@  $(LDFLAGS)

mpileup:
	$(CC) $(CFLAGS) -I $(HTS_DIR) mpileup.c -o $@  $(LDFLAGS)

modstate:
	$(CC) $(CFLAGS) -I $(HTS_DIR) modstate.c -o $@  $(LDFLAGS)

pileup_mod:
	$(CC) $(CFLAGS) -I $(HTS_DIR) pileup_mod.c -o $@  $(LDFLAGS)

flags_field:
	$(CC) $(CFLAGS) -I $(HTS_DIR) flags_htsopt_field.c -o $@  $(LDFLAGS)

split_t1:
	$(CC) $(CFLAGS) -I $(HTS_DIR) split_thread1.c -o $@  $(LDFLAGS)

split_t2:
	$(CC) $(CFLAGS) -I $(HTS_DIR) split_thread2.c -o $@  $(LDFLAGS)

clean:
	find . -name "*.o" | xargs rm -rf
	find . -name "*.dSYM" | xargs rm -rf
	rm $(PRGS)


