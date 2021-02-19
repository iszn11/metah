#pragma once

typedef struct Imetah
{
	void* (*init)(const void* config);
	void (*reset)(void* self);
	void (*deinit)(void* self);
	i32 (*iterate)(void* self);

	void (*fprint_detail_header)(void* self, FILE* file);
	void (*fprint_detail)(void* self, FILE* file);
	void (*fprint_overall_best)(void* self, FILE* file);
} Imetah;
