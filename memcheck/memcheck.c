#include <stdlib.h>
#include <stdio.h>
#include "htslib/vcf.h"
#include <assert.h>
main(int argc, char** argv)
{
  assert(argc >= 2 && "Requires 1 argument <vcf.gz> file name");
  //open file and read record
  htsFile* fptr = hts_open(argv[1], "r");
  assert(fptr);
  bcf_hdr_t* hdr = bcf_hdr_read(fptr);
  assert(hdr);
  int i = 0;
  bcf1_t* bcf_rec = 0;
  bcf1_t* duplicate = 0;
  int num_values = 0;
  int* read_end_value = (int*)malloc(4*sizeof(int));
  int original_end_point = -1;
  int modified_end_point = -1;
  int duplicate_end_point = -1;
  int status = -1;

  bcf_rec = bcf_init();
  bcf_read(fptr, hdr, bcf_rec);
  bcf_unpack(bcf_rec, BCF_UN_STR);
  bcf_update_id(hdr, bcf_rec, 0);
  //Create duplicate
  duplicate = bcf_dup(bcf_rec);
  //Set END for original record
  int new_end_point = 568;
  printf("Modifying END value of original record to %d\n",new_end_point);
  assert(bcf_update_info_int32(hdr, bcf_rec, "END", &new_end_point, 1) >= 0);
  /*//Re-read END for original record to see if it was correctly updated*/
  status = bcf_get_info_int32(hdr, bcf_rec, "END", &read_end_value,&num_values);
  assert(status == -3 || (status >=0 && num_values == 1));//either not present, or 1 element of correct type (int32)
  if(status >= 0)
    modified_end_point = read_end_value[0];
  else
    modified_end_point = bcf_rec->pos + 1;       //pos is 0 based, END is 1 based
  printf("Updated original record's end point to %d\n",modified_end_point);
  //Get END for duplicate record - NOTE: duplicate HAS NOT been modified
  bcf_unpack(duplicate, BCF_UN_INFO);
  status = bcf_get_info_int32(hdr, duplicate, "END", &read_end_value,&num_values);
  assert(status == -3 || (status >=0 && num_values == 1));//either not present, or 1 element of correct type (int32)
  if(status >= 0)
    duplicate_end_point = read_end_value[0];
  else
    duplicate_end_point = duplicate->pos + 1;       //pos is 0 based, END is 1 based
  printf("Re-reading unmodified duplicate record's end point : %d\n",duplicate_end_point);
  if(original_end_point != new_end_point && duplicate_end_point == new_end_point)
    printf("Eh? Why did the duplicate's end point change?\n");

  FILE* debug_fptr = stdout;
  kstring_t debug_string = { 0, 0, 0 };
  vcf_format(hdr,bcf_rec,&debug_string);
  fprintf(debug_fptr,"Modified original:: %s",debug_string.s);
  debug_string.l = 0;
  vcf_format(hdr,duplicate,&debug_string);
  fprintf(debug_fptr,"Un-modified duplicate:: %s",debug_string.s);
  fflush(debug_fptr);
  free(debug_string.s);
  debug_string.m = 0;

  bcf_destroy(bcf_rec);
  bcf_destroy(duplicate);

  bcf_hdr_destroy(hdr);
  hts_close(fptr);
  free(read_end_value);
  return 0;
}
