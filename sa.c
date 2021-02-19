#include "stdafx.h"
#include "sa.h"
#include "graph.h"

sa* sa_init(const sa_config* const config)
{
	ASSERT(config->initial_temp > 1.0 && config->neighbor_count > 0);

	const i16 gene_count = 2 * config->graph->size - 3;

	const i32 current_genes_offset = sizeof(sa);
	const i32 neighbor_genes_offset = current_genes_offset + gene_count * sizeof(i8);
	const i32 candidate_genes_offset = neighbor_genes_offset + gene_count * sizeof(i8);
	const i32 best_overall_genes_offset = candidate_genes_offset + gene_count * sizeof(i8);
	const i32 data_size = best_overall_genes_offset + gene_count * sizeof(i8);

	u8* mem = alloc(data_size);

	sa* const s = (sa*)mem;
	s->current_genes = (i8*)(mem + current_genes_offset);
	s->neighbor_genes = (i8*)(mem + neighbor_genes_offset);
	s->candidate_genes = (i8*)(mem + candidate_genes_offset);
	s->best_overall_genes = (i8*)(mem + best_overall_genes_offset);
	s->config = config;

	return s;
}

void sa_reset(sa* const s)
{
	const sa_config* const config = s->config;

	s->temp = config->initial_temp;

	s->config->init_fn(s->config->graph, s->current_genes, 1);
	s->current_fit = graph_fit_single(config->graph, s->current_genes);
	s->overall_best = s->current_fit;
	memcpy(s->best_overall_genes, s->current_genes, (2ULL * config->graph->size - 3) * sizeof(i8));

	s->iter = 0;
}

void sa_deinit(sa* const s)
{
	dealloc(s);
}

i32 sa_iterate(sa* const s)
{
	const sa_config* const config = s->config;
	const i32 neighbor_count      = config->neighbor_count;
	Imutate* const neighbor_fn    = config->neighbor_fn;
	i8* const current_genes       = s->current_genes;
	i8* const neighbor_genes      = s->neighbor_genes;
	i8* const candidate_genes     = s->candidate_genes;
	const i16 gene_count          = 2 * config->graph->size - 3;

	i32 candidate_fit = 0x7FFFFFFF;

	for (i32 i_neighbor = 0; i_neighbor < neighbor_count; ++i_neighbor)
	{
		mutate_copy(neighbor_fn, gene_count, current_genes, neighbor_genes);
		const i32 fit = graph_fit_single(config->graph, neighbor_genes);

		if (fit < candidate_fit)
		{
			memcpy(candidate_genes, neighbor_genes, gene_count * sizeof(i8));
			candidate_fit = fit;
		}
	}

	const float f = (float)exp(((double)s->current_fit - candidate_fit) / s->temp);
	if (candidate_fit < s->current_fit || random_chance(f))
	{
		s->current_fit = candidate_fit;
		s->current_genes = candidate_genes;
		s->candidate_genes = current_genes;
	}

	if (candidate_fit < s->overall_best)
	{
		s->overall_best = candidate_fit;
		memcpy(s->best_overall_genes, candidate_genes, gene_count * sizeof(i8));
	}

	s->temp = s->config->cooling_impl(s->config->cooling, s->temp, s->iter);
	s->iter += 1;

	return s->overall_best;
}

static void sa_Imetah_fprint_detail_header(sa* const s, FILE* file)
{
	(void)s;
	fprintf(file, "current\toverall_best\n");
}

static void sa_Imetah_fprint_detail(sa* const s, FILE* file)
{
	fprintf(file, "%d\t%d\t%d\n",
		s->iter,
		convert_units(s->current_fit),
		convert_units(s->overall_best)
	);
}

static void sa_Imetah_fprint_overall_best(sa* const s, FILE* file)
{
	const i16 gene_count = s->config->graph->size * 2 - 3;

	const i8* const best_genes = s->best_overall_genes;
	for (i16 i = 0; i < gene_count - 1; ++i)
	{
		fprintf(file, "%d ", best_genes[i]);
	}
	fprintf(file, "%d\n", best_genes[gene_count - 1]);
}

const Imetah sa_Imetah = {
	.init = sa_init,
	.reset = sa_reset,
	.deinit = sa_deinit,
	.iterate = sa_iterate,
	.fprint_detail_header = sa_Imetah_fprint_detail_header,
	.fprint_detail = sa_Imetah_fprint_detail,
	.fprint_overall_best = sa_Imetah_fprint_overall_best
};
