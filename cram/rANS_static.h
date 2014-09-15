#ifndef RANS_STATIC_H
#define RANS_STATIC_H

unsigned char *rans_compress(unsigned char *in, unsigned int in_size,
			     unsigned int *out_size, int order);
unsigned char *rans_uncompress(unsigned char *in, unsigned int in_size,
			       unsigned int *out_size, int order);


#endif /* RANS_STATIC_H */
