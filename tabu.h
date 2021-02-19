#pragma once
#include "operators.h"
#include "metah.h"

typedef struct graph graph;

typedef struct tabu_config
{
	const graph* graph;

	Iinit* init_fn;
	Imutate* neighbor_fn;
	Ifithash* fithash_fn;

	i32 tabu_capacity;
	i32 neighbor_count;
} tabu_config;

typedef struct tabu
{
	i8* current_genes;
	i8* neighbor_genes;
	i8* candidate_genes;
	i8* best_overall_genes;

	const tabu_config* config;

	double last_avg;

	u32* tabu_list;
	i32 tabu_length;
	i32 tabu_tail;

	i32 last_best;
	i32 last_worst;
	i32 overall_best;

	i32 iter;
} tabu;

tabu* tabu_init(const tabu_config* config);
void tabu_reset(tabu* t);
void tabu_deinit(tabu* t);

i32 tabu_iterate(tabu* t);

extern const Imetah tabu_Imetah;
