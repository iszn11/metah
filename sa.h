#pragma once
#include "operators.h"
#include "metah.h"

typedef struct graph graph;

typedef struct sa_config
{
	const graph* graph;

	Iinit* init_fn;
	Imutate* neighbor_fn;
	const void* cooling;
	Icooling* cooling_impl;

	double initial_temp;
	i32 neighbor_count;
} sa_config;

typedef struct sa
{
	i8* current_genes;
	i8* neighbor_genes;
	i8* candidate_genes;
	i8* best_overall_genes;

	const sa_config* config;

	double temp;

	i32 current_fit;
	i32 overall_best;

	i32 iter;
} sa;

sa* sa_init(const sa_config* config);
void sa_reset(sa* s);
void sa_deinit(sa* s);

i32 sa_iterate(sa* s);

extern const Imetah sa_Imetah;
