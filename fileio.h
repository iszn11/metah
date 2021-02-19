#pragma once

u8* file_read_alloc(const char* filename, void* (*alloc_fn)(i32));
