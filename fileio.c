#include "stdafx.h"

u8* file_read_alloc(const char* const filename, void* const (*alloc_fn)(i32))
{
	FILE* file;
	fopen_s(&file, filename, "rb");
	if (!file) abort();

	fseek(file, 0, SEEK_END);
	i32 file_size = ftell(file);
	rewind(file);

	if (file_size >= 0x7FFFFFFF) abort();
	u8* buf = alloc_fn(file_size + 1);

	if (fread(buf, 1, file_size, file) != file_size) abort();
	fclose(file);
	buf[file_size] = 0;

	return buf;
}
