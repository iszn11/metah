#include "stdafx.h"
#include "pop.h"
#include "graph.h"

pop* pop_init(const pop_config* const config)
{
	ASSERT(config->pop_size > 0);

	const i16 gene_count = 2 * config->graph->size - 3;
	const i16 pop_size = config->pop_size;

	const i32 fitness_offset = sizeof(pop);
	const i32 genes_offset = fitness_offset + pop_size * sizeof(i32);
	const i32 old_genes_offset = genes_offset + pop_size * gene_count * (i32)sizeof(i8);
	const i32 data_size = old_genes_offset + pop_size * gene_count * (i32)sizeof(i8);

	u8* mem = alloc(data_size);

	pop* const p = (pop*)mem;
	p->fitness = (i32*)(mem + fitness_offset);
	p->genes = (i8*)(mem + genes_offset);
	p->old_genes = (i8*)(mem + old_genes_offset);
	p->config = config;
	p->gene_count = gene_count;

	config->select_impl->select_data_init(config->select, p);

	return p;
}

void pop_reset(pop* const p)
{
	const pop_config* const config = p->config;

	config->init_fn(config->graph, p->genes, config->pop_size);
	if (config->force_canonical)
	{
		make_canonical_many(config->graph, p->genes, config->pop_size);
	}

	graph_fit_many(config->graph, p->genes, p->fitness, config->pop_size, &p->best_ind, &p->best_fitness, &p->worst_fitness);

	p->iter = 0;
}

void pop_deinit(pop* const p)
{
	const pop_config* const config = p->config;
	config->select_impl->select_data_deinit(config->select, p);
	dealloc(p);
}

i32 pop_iterate(pop* const p)
{
	const i16 gene_count             = p->gene_count;
	const pop_config* const config   = p->config;
	const float cross_chance         = config->cross_chance;
	const i16 pop_size               = config->pop_size;
	const void* const select         = config->select;
	const Iselect* const select_impl = config->select_impl;
	Icross* const cross_fn           = config->cross_fn;

	// swap genes
	i8* const tmp = p->genes;
	p->genes = p->old_genes;
	p->old_genes = tmp;

	// select
	select_impl->preprocess(select, p);

	i8* dest_genes = p->genes;
	for (i16 i = 0; i < pop_size; ++i, dest_genes += gene_count)
	{
		const i8* const parent1 = select_impl->select(select, p);

		if (random_chance(cross_chance))
		{
			const i8* const parent2 = select_impl->select(select, p);
			cross_fn(gene_count, parent1, parent2, dest_genes);
		}
		else
		{
			memcpy(dest_genes, parent1, gene_count * sizeof(i8));
		}
	}

	// mutate
	config->mutate_fn(p->gene_count, p->genes, pop_size, config->mutate_chance);

	// make canonical
	if (config->force_canonical)
	{
		make_canonical_many(config->graph, p->genes, pop_size);
	}

	// fitness
	graph_fit_many(config->graph, p->genes, p->fitness, pop_size, &p->best_ind, &p->best_fitness, &p->worst_fitness);

	p->iter += 1;
	return p->best_fitness;
}

static void pop_Imetah_fprint_detail_header(pop* const p, FILE* file)
{
	(void)p;
	fprintf(file, "iter\tbest\tworst\tavg\n");
}

static void pop_Imetah_fprint_detail(pop* const p, FILE* file)
{
	const i16 pop_size = p->config->pop_size;

	i32 sum = 0;
	for (i16 i = 0; i < pop_size; i++)
	{
		sum += p->fitness[i];
	}
	const double avg = (double)sum / pop_size;

	fprintf(file, "%d\t%d\t%d\t%d\n", p->iter, convert_units(p->best_fitness), convert_units(p->worst_fitness), convert_units(avg));
}

static void pop_Imetah_fprint_overall_best(pop* const p, FILE* file)
{
	const i16 gene_count = p->gene_count;

	const i8* const best_genes = &p->genes[p->best_ind * gene_count];
	for (i16 i = 0; i< gene_count - 1; ++i)
	{
		fprintf(file, "%d ", best_genes[i]);
	}
	fprintf(file, "%d\n", best_genes[gene_count - 1]);
}

const Imetah pop_Imetah = {
	.init = pop_init,
	.reset = pop_reset,
	.deinit = pop_deinit,
	.iterate = pop_iterate,
	.fprint_detail_header = pop_Imetah_fprint_detail_header,
	.fprint_detail = pop_Imetah_fprint_detail,
	.fprint_overall_best = pop_Imetah_fprint_overall_best
};
