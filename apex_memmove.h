/*
 * Apex memmove of T. Herselman.
 * The fastest memcpy/memmove on x86/x64 ... EVER,
 * written in C (according to the author).
 *
 * The original version has some fallback functions to use SSE/SSE2/SSE3,
 * I only keep the SSE4.2 version since most computers now support SSE4.2
 */


#ifndef _APEX_MEMMOVE_H_
#define _APEX_MEMMOVE_H_


extern void *(*apex_memcpy)(void *dst, const void *src, size_t size);
extern void *(*apex_memmove)(void *dst, const void *src, size_t size);


#endif
