#include "stdafx.h"
#include "tabu.h"
#include "graph.h"

static bool tabu_hash_exists(const tabu* t, u32 hash);
static void tabu_insert_hash(tabu* t, u32 hash);

tabu* tabu_init(const tabu_config* const config)
{
	ASSERT(config->tabu_capacity > 0 && config->neighbor_count > 0);

	const i16 gene_count = 2 * config->graph->size - 3;

	const i32 tabu_list_offset = sizeof(tabu);
	const i32 current_genes_offset = tabu_list_offset + config->tabu_capacity * sizeof(u32);
	const i32 neighbor_genes_offset = current_genes_offset + gene_count * sizeof(i8);
	const i32 candidate_genes_offset = neighbor_genes_offset + gene_count * sizeof(i8);
	const i32 best_overall_genes_offset = candidate_genes_offset + gene_count * sizeof(i8);
	const i32 data_size = best_overall_genes_offset + gene_count * sizeof(i8);

	u8* mem = alloc(data_size);

	tabu* const t = (tabu*)mem;
	t->current_genes = (i8*)(mem + current_genes_offset);
	t->neighbor_genes = (i8*)(mem + neighbor_genes_offset);
	t->candidate_genes = (i8*)(mem + candidate_genes_offset);
	t->best_overall_genes = (i8*)(mem + best_overall_genes_offset);
	t->config = config;
	t->tabu_list = (u32*)(mem + tabu_list_offset);

	return t;
}

void tabu_reset(tabu* const t)
{
	const tabu_config* const config = t->config;

	t->tabu_length = 0;
	t->tabu_tail = 0;

	config->init_fn(config->graph, t->current_genes, 1);
	t->overall_best = graph_fit_single(config->graph, t->current_genes);
	memcpy(t->best_overall_genes, t->current_genes, (2ULL * config->graph->size - 3) * sizeof(i8));

	t->iter = 0;
}

void tabu_deinit(tabu* const t)
{
	dealloc(t);
}

i32 tabu_iterate(tabu* t)
{
	const tabu_config* const config = t->config;
	const i32 neighbor_count = config->neighbor_count;
	Imutate* const neighbor_fn = config->neighbor_fn;
	Ifithash* const fithash_fn = config->fithash_fn;
	i8* const current_genes = t->current_genes;
	i8* const neighbor_genes = t->neighbor_genes;
	i8* const candidate_genes = t->candidate_genes;
	const i16 gene_count = 2 * config->graph->size - 3;

	double sum = 0.0;
	i32 last_best = 0x7FFFFFFF;
	i32 last_worst = 0;
	i32 candidate_fit = 0x7FFFFFFF;

	// NOTE initialize only to remove "potentially uninitialized" warning at tabu_insert_hash
	u32 candidate_hash = 0;

	for (i32 i_neighbor = 0; i_neighbor < neighbor_count; ++i_neighbor)
	{
		mutate_copy(neighbor_fn, gene_count, current_genes, neighbor_genes);
		i32 fit;
		u32 hash;
		fithash_fn(t, &fit, &hash);

		if (fit < candidate_fit && !tabu_hash_exists(t, hash))
		{
			memcpy(candidate_genes, neighbor_genes, gene_count * sizeof(i8));
			candidate_fit = fit;
			candidate_hash = hash;
		}

		sum += fit;
		if (fit < last_best) last_best = fit;
		if (fit > last_worst) last_worst = fit;
	}

	if (candidate_fit == 0x7FFFFFFF)
	{
		mutate_copy(neighbor_fn, gene_count, current_genes, candidate_genes);
		i32 fit;
		u32 hash;
		fithash_fn(t, &fit, &hash);
		candidate_fit = fit;
		candidate_hash = hash;
	}

	tabu_insert_hash(t, candidate_hash);

	t->current_genes = candidate_genes;
	t->candidate_genes = current_genes;
	t->last_avg = sum / neighbor_count;
	t->last_best = last_best;
	t->last_worst = last_worst;

	if (candidate_fit < t->overall_best)
	{
		t->overall_best = last_best;
		memcpy(t->best_overall_genes, candidate_genes, gene_count * sizeof(i8));
	}

	t->iter += 1;

	return t->overall_best;
}

static bool tabu_hash_exists(const tabu* const t, const u32 hash)
{
	const i32 tabu_length = t->tabu_length;
	const u32* const tabu_list = t->tabu_list;

	for (i32 i = 0; i < tabu_length; ++i)
	{
		const u32 h = tabu_list[i];
		if (hash == h)
		{
			return true;
		}
	}

	return false;
}

static void tabu_insert_hash(tabu* t, u32 hash)
{
	const i32 tabu_capacity = t->config->tabu_capacity;
	const i32 tabu_length = t->tabu_length;

	t->tabu_list[t->tabu_tail] = hash;
	t->tabu_length = min(tabu_capacity, tabu_length + 1);
	t->tabu_tail = (t->tabu_tail + 1) % tabu_capacity;
}

static void tabu_Imetah_fprint_detail_header(tabu* const t, FILE* file)
{
	(void)t;
	fprintf(file, "iter\tcurrent\tworst\tavg\toverall_best\n");
}

static void tabu_Imetah_fprint_detail(tabu* const t, FILE* file)
{
	fprintf(file, "%d\t%d\t%d\t%d\t%d\n",
		t->iter,
		convert_units(graph_fit_single(t->config->graph, t->current_genes)),
		convert_units(t->last_worst),
		convert_units(t->last_avg),
		convert_units(t->overall_best)
	);
}

static void tabu_Imetah_fprint_overall_best(tabu* const t, FILE* file)
{
	const i16 gene_count = t->config->graph->size * 2 - 3;

	const i8* const best_genes = t->best_overall_genes;
	for (i16 i = 0; i < gene_count - 1; ++i)
	{
		fprintf(file, "%d ", best_genes[i]);
	}
	fprintf(file, "%d\n", best_genes[gene_count - 1]);
}

const Imetah tabu_Imetah = {
	.init = tabu_init,
	.reset = tabu_reset,
	.deinit = tabu_deinit,
	.iterate = tabu_iterate,
	.fprint_detail_header = tabu_Imetah_fprint_detail_header,
	.fprint_detail = tabu_Imetah_fprint_detail,
	.fprint_overall_best = tabu_Imetah_fprint_overall_best
};
