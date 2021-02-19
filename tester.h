#pragma once

typedef struct Imetah Imetah;
typedef struct graph graph;

void test_metah(FILE* file, const void* config, const Imetah* metah_impl, i32 iterations, i32 repeats);
void test_random(FILE* file, const graph* graph, i32 count, i32 repeats);
void test_greedy(FILE* file, const graph* graph);
void test_detail(const void* config, const Imetah* metah_impl, i32 iterations, const char* filename);
