#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

static char* fuzzer_get_tmpfile(const uint8_t* data, size_t size) {
  char* filename_buffer = strdup("/tmp/generate_temporary_file.XXXXXX");
  if (!filename_buffer) {
    perror("Failed to allocate file name buffer.");
    abort();
  }
  const int file_descriptor = mkstemp(filename_buffer);
  if (file_descriptor < 0) {
    perror("Failed to make temporary file.");
    abort();
  }
  FILE* file = fdopen(file_descriptor, "wb");
  if (!file) {
    perror("Failed to open file descriptor.");
    close(file_descriptor);
    abort();
  }
  const size_t bytes_written = fwrite(data, sizeof(uint8_t), size, file);
  if (bytes_written < size) {
    close(file_descriptor);
    fprintf(stderr, "Failed to write all bytes to file (%zu out of %zu)",
            bytes_written, size);
    abort();
  }
  fclose(file);
  return filename_buffer;
}

static void fuzzer_release_tmpfile(char* filename) {
  if (unlink(filename) != 0) {
    perror("WARNING: Failed to delete temporary file.");
  }
  free(filename);
}


// Duplicated from: htsfile.c
static htsFile *dup_stdout(const char *mode) {
  int fd = dup(STDOUT_FILENO);
  hFILE *hfp = (fd >= 0) ? hdopen(fd, mode) : NULL;
  return hfp ? hts_hopen(hfp, "-", mode) : NULL;
}

static int view_sam(htsFile *in) {
  if (!in) return 0;
  samFile *out = dup_stdout("w");
  bam_hdr_t *hdr = sam_hdr_read(in);

  (void)sam_hdr_write(out, hdr);
  bam1_t *b = bam_init1();
  while (sam_read1(in, hdr, b) >= 0) (void)sam_write1(out, hdr, b);
  bam_destroy1(b);

  bam_hdr_destroy(hdr);
  hts_close(out);
  return 1;
}

static int view_vcf(htsFile *in) {
  if (!in) return 0;
  vcfFile *out = dup_stdout("w");
  bcf_hdr_t *hdr = bcf_hdr_read(in);

  bcf_hdr_write(out, hdr);
  bcf1_t *rec = bcf_init();
  while (bcf_read(in, hdr, rec) >= 0) bcf_write(out, hdr, rec);
  bcf_destroy(rec);

  bcf_hdr_destroy(hdr);
  hts_close(out);
  return 1;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  char* temporary_file = fuzzer_get_tmpfile(data, size);
  htsFile *ht_file = hts_open(temporary_file, /*mode=*/"r");
  if (ht_file == NULL) {
    fuzzer_release_tmpfile(temporary_file);
    return 0;
  }
  switch (ht_file->format.category) {
    case sequence_data:
      view_sam(ht_file);
      break;
    case variant_data:
      view_vcf(ht_file);
      break;
    default:
      break;
  }
  hts_close(ht_file);
  fuzzer_release_tmpfile(temporary_file);
  return 0;
}
