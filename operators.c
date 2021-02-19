#include "stdafx.h"
#include "graph.h"
#include "operators.h"
#include "pop.h"
#include "tabu.h"

// --- INITIALIZATION ----------------------------------------------------------

void random_Iinit(const graph* const graph, i8* genes, const i32 count)
{
	const i16 city_count = graph->size - 1;
	const i16 gene_count = city_count * 2 - 1;

	for (i32 i = 0; i < count; ++i, genes += gene_count)
	{
		i16 g;
		for (g = 0; g < city_count; ++g) genes[g] = (i8)(g + 1);
		for (; g < gene_count; ++g) genes[g] = 0;
		partial_shuffle_i8(genes, gene_count - 1, gene_count);
	}
}

void greedy_Iinit(const graph* const graph, i8* genes, const i32 count)
{
	const i16 city_count = graph->size - 1;
	const i16 gene_count = city_count * 2 - 1;
	const i8 max_capacity = graph->capacity;

	for (i16 i_pop = 0; i_pop < count; ++i_pop, genes += gene_count)
	{
		i16 g;
		for (g = 0; g < city_count; ++g) genes[g] = (i8)(g + 1);
		for (; g < gene_count; ++g) genes[g] = 0;

		i8 capacity_left = max_capacity;
		i8 from = 0;

		i16 gene_ind = 0;
		i16 past_city_end_ind = city_count;

		for (i8 i_city = 0; i_city < city_count; ++i_city, ++gene_ind)
		{
			if (from == 0)
			{
				if (i_city < city_count - 1)
				{
					partial_shuffle_i8(&genes[gene_ind], 1, past_city_end_ind - gene_ind);
				}
			}
			else
			{
				i16 shortest_to_ind = gene_ind;
				i16 shortest_length = graph_get_weight(graph, from, genes[shortest_to_ind]);
				for (i16 edge_ind = gene_ind + 1; edge_ind < past_city_end_ind; ++edge_ind)
				{
					const i16 edge_length = graph_get_weight(graph, from, genes[edge_ind]);
					if (edge_length < shortest_length)
					{
						shortest_to_ind = edge_ind;
						shortest_length = edge_length;
					}
				}

				const i8 tmp = genes[gene_ind];
				genes[gene_ind] = genes[shortest_to_ind];
				genes[shortest_to_ind] = tmp;
			}

			const i8 to = genes[gene_ind];
			const i8 demand = graph->demands[to];

			capacity_left -= demand;
			if (capacity_left < 0)
			{
				ASSERT(past_city_end_ind < gene_count);
				genes[past_city_end_ind] = genes[gene_ind];
				genes[gene_ind] = 0;

				past_city_end_ind += 1;
				i_city -= 1;

				capacity_left = max_capacity;
				from = 0;
			}
			else
			{
				from = to;
			}
		}
	}
}

// --- SELECTION ---------------------------------------------------------------

static i8* tourney_Iselect_select_distinct(const void* const self, const pop* const p)
{
	const tourney* const tourney   = self;
	const i32* const fitness       = p->fitness;
	i8* const old_genes            = p->old_genes;
	const i16 gene_count           = p->gene_count;
	const pop_config* const config = p->config;
	const i16 pop_size             = config->pop_size;
	const float f_size             = tourney->size;
	ASSERT(f_size >= 0.0f && f_size <= 1.0f);
	const i16 i_size               = (i16)max(pop_size * f_size, 1);

	i16* buf = alloc_tmp(pop_size * sizeof(i16));
	for (i16 i = 0; i < pop_size; ++i) buf[i] = i;
	partial_shuffle_i16(buf, i_size, pop_size);

	i32 best = 0x7FFFFFFF;
	i8* ret = NULL;

	for (i16 i = 0; i < i_size; ++i)
	{
		const i32 f = fitness[buf[i]];
		if (f < best)
		{
			best = f;
			ret = &old_genes[buf[i] * gene_count];
		}
	}

	return ret;
}

static i8* tourney_Iselect_select(const void* const self, const pop* const p)
{
	const tourney* const tourney   = self;
	const i32* const fitness       = p->fitness;
	i8* const old_genes            = p->old_genes;
	const i16 gene_count           = p->gene_count;
	const pop_config* const config = p->config;
	const i16 pop_size             = config->pop_size;
	const float f_size             = tourney->size;
	ASSERT(f_size >= 0.0f && f_size <= 1.0f);
	const i16 i_size               = (i16)max(pop_size * f_size, 1);

	i32 best = 0x7FFFFFFF;
	i8* ret = NULL;

	for (i16 i = 0; i < i_size; ++i)
	{
		i16 ind = (i16)random_range(0, pop_size);
		const i32 f = fitness[ind];
		if (f < best)
		{
			best = f;
			ret = &old_genes[ind * gene_count];
		}
	}

	return ret;
}

static void tourney_Iselect_preprocess(const void* const self, const pop* const src) { (void)self; (void)src; }
static void tourney_Iselect_select_data_init(const void* const self, pop* const p) { (void)self; (void)p; }
static void tourney_Iselect_select_data_deinit(const void* const self, pop* const p) { (void)self; (void)p; }

const Iselect tourney_Iselect = {
	.select = tourney_Iselect_select,
	.preprocess = tourney_Iselect_preprocess,
	.select_data_init = tourney_Iselect_select_data_init,
	.select_data_deinit = tourney_Iselect_select_data_deinit,
};

const Iselect tourney_Iselect_distinct = {
	.select = tourney_Iselect_select_distinct,
	.preprocess = tourney_Iselect_preprocess,
	.select_data_init = tourney_Iselect_select_data_init,
	.select_data_deinit = tourney_Iselect_select_data_deinit,
};

static i8* roulette_Iselect_select(const void* const self, const pop* const p)
{
	(void)self;
	const float* weight = p->select_data;

	float w = random_float();
	i16 ind = 0;
	while (w >= 0.0f)
	{
		w -= *weight;
		weight += 1;
		ind += 1;
	}

	return &p->old_genes[(ind - 1) * p->gene_count];
}

static void roulette_Iselect_preprocess(const void* const self, const pop* const p)
{
	const roulette* const roulette = self;
	const float k                  = roulette->k;
	ASSERT(k > 0.0f);
	const i16 pop_size             = p->config->pop_size;
	const i32* fitness             = p->fitness;
	float* weight                  = p->select_data;
	const i32 min                  = p->best_fitness;
	const i32 max                  = p->worst_fitness;

	const i32 d = max - min;
	double sum = 0.0;
	for (i16 i = 0; i < pop_size; ++i, ++fitness, ++weight)
	{
		const i32 f = *fitness;
		const double w = pow(((double)max - f) / d, k);
		*weight = (float)w;
		sum += w;
	}

	weight = p->select_data;
	for (i16 i = 0; i < pop_size - 1; ++i, ++weight)
	{
		*weight = (float)(*weight / sum);
	}
	weight[pop_size - 1] = 1.0f;
}

static void roulette_Iselect_select_data_init(const void* const self, pop* const p)
{
	(void)self;
	p->select_data = alloc(p->config->pop_size * sizeof(float));
}

static void roulette_Iselect_select_data_deinit(const void* const self, pop* const p)
{
	(void)self;
	dealloc(p->select_data);
	p->select_data = NULL;
}

const Iselect roulette_Iselect = {
	.select = roulette_Iselect_select,
	.preprocess = roulette_Iselect_preprocess,
	.select_data_init = roulette_Iselect_select_data_init,
	.select_data_deinit = roulette_Iselect_select_data_deinit,
};

// --- CROSS -------------------------------------------------------------------

void ox_Icross(const i16 gene_count, const i8* const parent1, const i8* const parent2, i8* const child) // TODO ignore 0
{
	i16 start = (i16)random_range(0, gene_count);
	i16 end = (i16)random_range(0, gene_count);

	if (start > end)
	{
		const i16 tmp = start;
		start = end;
		end = tmp;
	}

	memcpy(&child[start], &parent1[start], ((i64)end - start + 1) * sizeof(i8));

	bool* used = alloc_tmp(gene_count * sizeof(bool));
	memset(used, 0, gene_count * sizeof(bool));
	for (i16 i = start; i <= end; ++i)
	{
		used[parent1[i]] = true;
	}

	i16 zero_count = 0;
	const i16 max_zero_count = gene_count / 2;
	for (i16 i = start; i <= end; ++i) if (child[i] == 0) zero_count += 1;

	i16 p_ind = 0;
	i16 c_ind = 0;

	for (; p_ind < gene_count; ++p_ind)
	{
		// skip past copied part from parent1
		if (c_ind >= start && c_ind <= end) c_ind = end + 1;

		const i8 g = parent2[p_ind];

		if (g == 0)
		{
			if (zero_count >= max_zero_count) continue;
			zero_count += 1;
		}
		else
		{
			if (used[g]) continue;
		}

		child[c_ind] = g;
		c_ind += 1;
	}
}

void pmx_Icross(i16 gene_count, const i8* parent1, const i8* parent2, i8* child)
{
	const i8 city_count = (i8)(gene_count + 1) / 2;

	i16 start = (i16)random_range(0, gene_count);
	i16 end = (i16)random_range(0, gene_count);

	if (start > end)
	{
		const i16 tmp = start;
		start = end;
		end = tmp;
	}

#ifdef VERBOSE
	printf("PMX CROSS\n");
	printf("P1: ");
	for (i16 i = 0; i < gene_count; ++i)
	{
		if (i == start) printf("| ");
		printf("%2d ", parent1[i]);
		if (i == end) printf("| ");
	}
	printf("\nP2: ");
	for (i16 i = 0; i < gene_count; ++i)
	{
		if (i == start) printf("| ");
		printf("%2d ", parent2[i]);
		if (i == end) printf("| ");
	}
#endif

	memcpy(&child[start], &parent1[start], ((i64)end - start + 1) * sizeof(i8));

	/* 0..city_count-1 maps cities 1..city_count
	 * city_count..gene_count maps zeros 1..city_count-1
	 */
	i16* const map = alloc_tmp(gene_count * sizeof(i16));
	for (i16 i = 0; i < gene_count; ++i) map[i] = i;
	{
		i16 i;
		i16 z1 = 0;
		i16 z2 = 0;
		for (i = 0; i < start; ++i)
		{
			const i8 g1 = parent1[i];
			const i8 g2 = parent2[i];

			if (g1 == 0) z1 += 1;
			if (g2 == 0) z2 += 1;
		}
		for (;i <= end; ++i)
		{
			const i8 g1 = parent1[i];
			const i8 g2 = parent2[i];

			if (g1 == 0) z1 += 1;
			if (g2 == 0) z2 += 1;

			const i16 map_target = (g2 == 0 ? city_count + z2 : g2) - 1;
			map[(g1 == 0 ? city_count + z1 : g1) - 1] = map_target;
		}
	}
#ifdef VERBOSE
	printf("\nMP: ");
	for (i16 i = 0; i < gene_count; ++i)
	{
		if (map[i] == i) continue;

		if (i >= city_count) printf("0.%d->", i - city_count + 1);
		else printf("%d->", i + 1);

		if (map[i] >= city_count) printf("0.%d ", map[i] - city_count + 1);
		else printf("%d ", map[i] + 1);
	}
#endif

	for (i16 i = 0, z = 0; i < gene_count; ++i)
	{
		const i8 g = parent2[i];
		if (g == 0) z += 1;

		// skip past copied part from parent1
		if (i >= start && i <= end) continue;

		i16 map_ind = (g == 0 ? city_count + z : g) - 1;
		while (map_ind != map[map_ind]) map_ind = map[map_ind];

		child[i] = map_ind >= city_count ? 0 : (i8)(map_ind + 1);
	}
#ifdef VERBOSE
	printf("\nCH: ");
	for (i16 i = 0; i < gene_count; ++i)
	{
		if (i == start) printf("| ");
		printf("%2d ", child[i]);
		if (i == end) printf("| ");
	}
	printf("\n\n");
#endif
}

// --- MUTATION ----------------------------------------------------------------

void swap_single_Imutate(const i16 gene_count, i8* genes, const i32 count, const float chance)
{
	ASSERT(chance >= 0.0f && chance <= 1.0f && count > 0);

	for (i16 i = 0; i < count; ++i, genes += gene_count)
	{
		if (!random_chance(chance))
		{
			continue;
		}

		const i16 g1 = (i16)random_range(0, gene_count);
		const i16 g2 = (i16)random_range(0, gene_count);

		const i8 tmp = genes[g1];
		genes[g1] = genes[g2];
		genes[g2] = tmp;
	}
}

void swap_many_Imutate(const i16 gene_count, i8* genes, const i32 count, const float chance)
{
	ASSERT(chance >= 0.0f && chance <= 1.0f && count > 0);

	for (i16 i = 0; i < count; ++i, genes += gene_count)
	{
		for (i16 j = 0; j < gene_count; ++j)
		{
			if (!random_chance(chance)) continue;

			const i16 other = (i16)random_range(0, gene_count);

			const i8 tmp = genes[j];
			genes[j] = genes[other];
			genes[other] = tmp;
		}
	}
}

void invert_Imutate(const i16 gene_count, i8* genes, const i32 count, const float chance)
{
	ASSERT(chance >= 0.0f && chance <= 1.0f);

	for (i16 i = 0; i < count; ++i, genes += gene_count)
	{
		if (!random_chance(chance))
		{
			continue;
		}

		i16 start = (i16)random_range(0, gene_count);
		i16 end = (i16)random_range(0, gene_count);
		if (start > end)
		{
			const i16 tmp = start;
			start = end;
			end = tmp;
		}

		for (i16 l = start, r = end; l < r; ++l, --r)
		{
			const i8 tmp = genes[l];
			genes[l] = genes[r];
			genes[r] = tmp;
		}
	}
}

void mutate_copy(Imutate* const op, const i16 gene_count, const i8* const src, i8* const dst)
{
	memcpy(dst, src, gene_count * sizeof(i8));
	op(gene_count, dst, 1, 1.0f);
}

// --- COOLING -----------------------------------------------------------------

double cool_Icooling_subtract(const void* self, double temp, i32 iter)
{
	(void)iter;
	const cool* const c = (const cool*)self;
	return fmax(temp - c->k, c->t_min);
}

double cool_Icooling_exp(const void* self, double temp, i32 iter)
{
	(void)iter;
	const cool* const c = (const cool*)self;
	const double t_min = c->t_min;
	return (temp - t_min) * c->k + t_min;
}

// --- FITHASH -----------------------------------------------------------------

void simple_Ifithash(const tabu* const t, i32* const fit, u32* const hash)
{
	const graph* const graph = t->config->graph;
	const i8 graph_size = graph->size;
#ifndef NDEBUG
	const i16 gene_count = 2 * graph_size - 3;
#endif
	const i16 city_count = graph_size - 1;
	const i8 max_capacity = graph->capacity;
	const i8* const genes = t->neighbor_genes;

	i32 f = 0;
	i16 cities = 0;

	i8 capacity_left = max_capacity;
	i8 from = 0;

	u32 h = 17;

	for (i16 j = 0; true; ++j)
	{
		ASSERT(j < gene_count);
		const i8 to = genes[j];
		const i8 demand = graph->demands[to];

		capacity_left -= demand;

		if (capacity_left < 0)
		{
			// NOTE to != 0, because we wouldn't ran out of capacity in that case
			// NOTE from != 0, because it would lead to an unsolvable case
			f += graph_get_weight(graph, from, 0);
			h = h * 23;
			f += graph_get_weight(graph, 0, to);
			h = h * 23 + to;
			capacity_left = max_capacity - demand;
			ASSERT(capacity_left >= 0);
		}
		else
		{
			f += graph_get_weight(graph, from, to);
			if (from != 0) h = h * 23 + to;
		}

		from = to;

		if (to != 0) cities += 1;
		if (cities >= city_count) break;
		if (to == 0) capacity_left = max_capacity;
	}
	f += graph_get_weight(graph, from, 0);

	*fit = f;
	*hash = h;
}

struct tabu_hash_journey
{
	i16 low_ind;
	i16 high_ind;
	i8 start_city;
	i8 dir;
};

static void sort_journeys(struct tabu_hash_journey* const tab, const i32 jc)
{
	for (i32 i = 1; i < jc; ++i)
	{
		struct tabu_hash_journey x = tab[i];
		i32 j = i - 1;
		while (j >= 0 && tab[j].start_city > x.start_city)
		{
			tab[j + 1] = tab[j];
			j -= 1;
		}
		tab[j + 1] = x;
	}
}

void complex_Ifithash(const tabu* const t, i32* const fit, u32* const hash)
{
	const graph* const graph = t->config->graph;
	const i8 graph_size = graph->size;
#ifndef NDEBUG
	const i16 gene_count = 2 * graph_size - 3;
#endif
	const i16 city_count = graph_size - 1;
	const i8 max_capacity = graph->capacity;
	const i8* const genes = t->neighbor_genes;

	i32 f = 0;
	i16 cities = 0;

	i8 capacity_left = max_capacity;
	i8 from = 0;

	struct tabu_hash_journey* j = alloc_tmp(city_count * sizeof(struct tabu_hash_journey));
	i16 jc = 0;
	bool on_journey = false;

	i16 i;
	for (i = 0; true; ++i)
	{
		ASSERT(i < gene_count);
		const i8 to = genes[i];
		const i8 demand = graph->demands[to];

		capacity_left -= demand;

		if (capacity_left < 0)
		{
			f += graph_get_weight(graph, from, 0);
			f += graph_get_weight(graph, 0, to);
			capacity_left = max_capacity - demand;
			ASSERT(capacity_left >= 0);

			ASSERT(on_journey);
			j[jc].high_ind = i - 1;

			const i8 low_city = genes[j[jc].low_ind];
			const i8 high_city = genes[j[jc].high_ind];
			ASSERT(low_city != 0 && high_city != 0);
			if (low_city < high_city)
			{
				j[jc].dir = 1;
				j[jc].start_city = low_city;
			}
			else
			{
				j[jc].dir = -1;
				j[jc].start_city = high_city;
			}

			on_journey = false;
			jc += 1;
		}
		else
		{
			f += graph_get_weight(graph, from, to);
		}

		if (to == 0)
		{
			capacity_left = max_capacity;
			if (from != 0)
			{
				ASSERT(on_journey);
				j[jc].high_ind = i - 1;

				const i8 low_city = genes[j[jc].low_ind];
				const i8 high_city = genes[j[jc].high_ind];
				ASSERT(low_city != 0 && high_city != 0);
				if (low_city < high_city)
				{
					j[jc].dir = 1;
					j[jc].start_city = low_city;
				}
				else
				{
					j[jc].dir = -1;
					j[jc].start_city = high_city;
				}

				on_journey = false;
				jc += 1;
			}
		}
		else
		{
			cities += 1;
			if (!on_journey)
			{
				j[jc].low_ind = i;
				on_journey = true;
			}
		}

		from = to;

		if (cities >= city_count) break;
	}
	f += graph_get_weight(graph, from, 0);
	ASSERT(on_journey);
	j[jc].high_ind = i;

	const i8 low_city = genes[j[jc].low_ind];
	const i8 high_city = genes[j[jc].high_ind];
	ASSERT(low_city != 0 && high_city != 0);
	if (low_city < high_city)
	{
		j[jc].dir = 1;
		j[jc].start_city = low_city;
	}
	else
	{
		j[jc].dir = -1;
		j[jc].start_city = high_city;
	}
	jc += 1;

	sort_journeys(j, jc);
	u32 h = 17;
	for (i8 i_journey = 0; i_journey < jc; ++i_journey)
	{
		const struct tabu_hash_journey* const jj = &j[i_journey];
		const i8 dir = jj->dir;

		if (dir == -1)
		{
			for (i16 ind = jj->high_ind; ind >= jj->low_ind; --ind)
			{
				const i8 city = genes[ind];
				ASSERT(city != 0);
				h = h * 23 + city;
			}
		}
		else
		{
			for (i16 ind = jj->low_ind; ind <= jj->high_ind; ++ind)
			{
				const i8 city = genes[ind];
				ASSERT(city != 0);
				h = h * 23 + city;
			}
		}

		h *= 23;
	}

	*fit = f;
	*hash = h;
}

void make_canonical_many(const graph* const graph, i8* genes, const i32 count)
{
	const i8 graph_size = graph->size;
	const i16 gene_count = 2 * graph_size - 3;
	for (i32 i = 0; i < count; ++i, genes += gene_count)
	{
		make_canonical_single(graph, genes);
	}
}

void make_canonical_single(const graph* const graph, i8* const genes)
{
	const i8 graph_size = graph->size;
	const i16 gene_count = 2 * graph_size - 3;
	const i16 city_count = graph_size - 1;
	const i8 max_capacity = graph->capacity;

	i16 cities = 0;

	i8 capacity_left = max_capacity;
	i8 from = 0;

	u8* const tmp = alloc_tmp(city_count * sizeof(struct tabu_hash_journey) + gene_count * sizeof(i8));
	const i8* const orig_genes = (const i8*)(tmp + city_count * sizeof(struct tabu_hash_journey));
	memcpy(tmp + city_count * sizeof(struct tabu_hash_journey), genes, gene_count * sizeof(i8));
	struct tabu_hash_journey* j = (struct tabu_hash_journey*)tmp;

#ifdef VERBOSE
	printf("\nGN: ");
	for (i16 i = 0; i < gene_count; ++i)
	{
		printf("%2d ", genes[i]);
	}
	printf("\nOR: ");
	for (i16 i = 0; i < gene_count; ++i)
	{
		printf("%2d ", orig_genes[i]);
	}
	printf("\n");
#endif

	i16 jc = 0;
	bool on_journey = false;

	i16 i;
	for (i = 0; true; ++i)
	{
		ASSERT(i < gene_count);
		const i8 to = orig_genes[i];
		const i8 demand = graph->demands[to];

		capacity_left -= demand;

		if (capacity_left < 0)
		{
			capacity_left = max_capacity - demand;
			ASSERT(capacity_left >= 0);

			ASSERT(on_journey);
			j[jc].high_ind = i - 1;

			const i8 low_city = orig_genes[j[jc].low_ind];
			const i8 high_city = orig_genes[j[jc].high_ind];
			ASSERT(low_city != 0 && high_city != 0);
			if (low_city < high_city)
			{
				j[jc].dir = 1;
				j[jc].start_city = low_city;
			}
			else
			{
				j[jc].dir = -1;
				j[jc].start_city = high_city;
			}

			on_journey = false;
			jc += 1;
		}

		if (to == 0)
		{
			capacity_left = max_capacity;
			if (from != 0)
			{
				ASSERT(on_journey);
				j[jc].high_ind = i - 1;

				const i8 low_city = orig_genes[j[jc].low_ind];
				const i8 high_city = orig_genes[j[jc].high_ind];
				ASSERT(low_city != 0 && high_city != 0);
				if (low_city < high_city)
				{
					j[jc].dir = 1;
					j[jc].start_city = low_city;
				}
				else
				{
					j[jc].dir = -1;
					j[jc].start_city = high_city;
				}

				on_journey = false;
				jc += 1;
			}
		}
		else
		{
			cities += 1;
			if (!on_journey)
			{
				j[jc].low_ind = i;
				on_journey = true;
			}
		}

		from = to;

		if (cities >= city_count) break;
	}

	ASSERT(on_journey);
	j[jc].high_ind = i;

	const i8 low_city = orig_genes[j[jc].low_ind];
	const i8 high_city = orig_genes[j[jc].high_ind];
	ASSERT(low_city != 0 && high_city != 0);
	if (low_city < high_city)
	{
		j[jc].dir = 1;
		j[jc].start_city = low_city;
	}
	else
	{
		j[jc].dir = -1;
		j[jc].start_city = high_city;
	}
	jc += 1;

	sort_journeys(j, jc);
	memset(genes, 0, gene_count * sizeof(i8));

	i8* g = genes;

	for (i8 i_journey = 0; i_journey < jc; ++i_journey)
	{
		const struct tabu_hash_journey* const jj = &j[i_journey];
		const i8 dir = jj->dir;

		if (dir == -1)
		{
			for (i16 ind = jj->high_ind; ind >= jj->low_ind; --ind)
			{
				const i8 city = orig_genes[ind];
				ASSERT(city != 0);
				*g++ = city;
			}
		}
		else
		{
			for (i16 ind = jj->low_ind; ind <= jj->high_ind; ++ind)
			{
				const i8 city = orig_genes[ind];
				ASSERT(city != 0);
				*g++ = city;
			}
		}

		if (i_journey - 1 < jc) *g++ = 0;
	}

#ifdef VERBOSE
	printf("\nGN: ");
	for (i16 ii = 0; ii < gene_count; ++ii)
	{
		printf("%2d ", genes[ii]);
	}
	printf("\n\n");
#endif
}
