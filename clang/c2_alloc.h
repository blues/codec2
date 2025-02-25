/* debug_alloc.h
 *
 * Some macros which can report on malloc results.
 *
 * Enable with "-D DEBUG_ALLOC"
 */

#ifndef DEBUG_ALLOC_H
#define DEBUG_ALLOC_H

#include <stdio.h>

// Debug calls

#ifdef CORTEX_M4
extern char *__heap_end;
register char *sp asm("sp");
#endif

extern void *codec2_malloc(size_t size);
extern void *codec2_calloc(size_t nmemb, size_t size);
extern void codec2_free(void *ptr);

#define MALLOC(size) codec2_malloc(size)
#define CALLOC(nmemb, size) codec2_calloc(nmemb, size)
#define FREE(ptr) codec2_free(ptr)

#endif  // DEBUG_ALLOC_H
