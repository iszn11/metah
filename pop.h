#pragma once
#include "operators.h"
#include "metah.h"

typedef struct graph graph;

typedef struct pop_config
{
	const graph* graph;

	Iinit* init_fn;
	const void* select;
	const Iselect* select_impl;
	Icross* cross_fn;
	Imutate* mutate_fn;

	float cross_chance;
	float mutate_chance;
	bool force_canonical;

	i16 pop_size;
} pop_config;

typedef struct pop
{
	i32* fitness;
	i8* genes;
	i8* old_genes;
	void* select_data;

	const pop_config* config;

	i32 best_fitness;
	i32 worst_fitness;
	i32 best_ind;
	i16 gene_count;

	i32 iter;
} pop;

pop* pop_init(const pop_config* config);
void pop_reset(pop* p);
void pop_deinit(pop* p);

i32 pop_iterate(pop* p);

extern const Imetah pop_Imetah;
