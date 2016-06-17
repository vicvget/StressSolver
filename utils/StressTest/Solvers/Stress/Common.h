#pragma once
#ifdef GNU_CPP
#define aligned_alloc(size, alignment) _mm_malloc(size, alignment)
#define aligned_free(ptr)	_mm_free(ptr)
#else
#define aligned_alloc(size,alignment) _aligned_malloc(size,alignment)
#define aligned_free(ptr)	_aligned_free(ptr)
#endif
	