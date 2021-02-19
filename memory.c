#include "stdafx.h"

void* alloc(const i32 size)
{
	ASSERT(size > 0);

	void* ptr = malloc(size);
	if (ptr == NULL) abort();

	return ptr;
}

void dealloc(void* const ptr)
{
	free(ptr);
}

static u8 tmp_buffer[1048576];

void* alloc_tmp(const i32 size)
{
	ASSERT(size > 0);

	if (size > sizeof(tmp_buffer)) abort();

	return tmp_buffer;
}
