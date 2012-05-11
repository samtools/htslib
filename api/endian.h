#ifndef ENDIAN_H
#define ENDIAN_H

#include <stdint.h>

static inline int ed_is_big()
{
	long one= 1;
	return !(*((char *)(&one)));
}
static inline uint16_t ed_swap_2(uint16_t v)
{
	return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *ed_swap_2p(void *x)
{
	*(uint16_t*)x = ed_swap_2(*(uint16_t*)x);
	return x;
}
static inline uint32_t ed_swap_4(uint32_t v)
{
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *ed_swap_4p(void *x)
{
	*(uint32_t*)x = ed_swap_4(*(uint32_t*)x);
	return x;
}
static inline uint64_t ed_swap_8(uint64_t v)
{
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *ed_swap_8p(void *x)
{
	*(uint64_t*)x = ed_swap_8(*(uint64_t*)x);
	return x;
}

#endif
