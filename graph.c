#include "stdafx.h"
#include "graph.h"
#include "pop.h"

static const char NAME_PREFIX[] = "NAME";
static const char SIZE_PREFIX[] = "DIMENSION";
static const char CAPACITY_PREFIX[] = "CAPACITY";
static const char NODE_SECTION[] = "NODE_COORD_SECTION";
static const char DEMAND_SECTION[] = "DEMAND_SECTION";

static const double SCALING = 141.42135623730950;

static graph* graph_init(i8 size, i32 name_length);

i32 convert_units(const double internal_dist)
{
	ASSERT(internal_dist >= 0.0);
	return lround(internal_dist * SCALING / 0x7FFF);
}

void graph_fit_many(
	const graph* const graph,
	const i8* genes,
	i32* fitness,
	const i32 count,
	i32* const best_ind,
	i32* const best_fit,
	i32* const worst_fit)
{
	const i8 graph_size = graph->size;
	const i16 gene_count = 2 * graph_size - 3;
	const i16 city_count = graph_size - 1;
	const i8 max_capacity = graph->capacity;

	i32 best_ind_value = -1;
	i32 best_fitness = 0x7FFFFFFF;
	i32 worst_fitness = 0;

	for (i32 i = 0; i < count; ++i, ++fitness, genes += gene_count)
	{
		i32 f = 0;
		i16 cities = 0;

		i8 capacity_left = max_capacity;
		i8 from = 0;

		for (i16 j = 0; true; ++j)
		{
			ASSERT(j < gene_count);
			const i8 to = genes[j];
			const i8 demand = graph->demands[to];

			capacity_left -= demand;

			if (capacity_left < 0)
			{
				f += graph_get_weight(graph, from, 0);
				f += graph_get_weight(graph, 0, to);
				capacity_left = max_capacity - demand;
				ASSERT(capacity_left >= 0);
			}
			else
			{
				f += graph_get_weight(graph, from, to);
			}

			from = to;

			if (to != 0) cities += 1;
			if (cities >= city_count) break;
			if (to == 0) capacity_left = max_capacity;
		}
		f += graph_get_weight(graph, from, 0);

		*fitness = f;

		if (f < best_fitness)
		{
			best_fitness = f;
			best_ind_value = i;
		}
		if (f > worst_fitness)
		{
			worst_fitness = f;
		}
	}

	if (best_ind) *best_ind = best_ind_value;
	if (best_fit) *best_fit = best_fitness;
	if (worst_fit) *worst_fit = worst_fitness;
}

i32 graph_fit_single(const graph* const graph, const i8* genes)
{
	const i16 city_count = graph->size - 1;
	const i8 max_capacity = graph->capacity;
	const i8* const demands = graph->demands;
	const i16* const weights = graph->weights;

	i32 f = 0;
	i16 cities = 0;

	i8 capacity_left = max_capacity;
	i8 from = 0;

	while (true)
	{
		const i8 to = *genes;
		const i8 demand = demands[to];

		if (capacity_left < demand)
		{
			f += weights[from];
			f += weights[to];
			capacity_left = max_capacity - demand;
			ASSERT(capacity_left >= 0);
		}
		else
		{
			f += graph_get_weight(graph, from, to);
			capacity_left -= demand;
		}

		from = to;

		if (to != 0) cities += 1;
		if (cities >= city_count) break;
		if (to == 0) capacity_left = max_capacity;

		++genes;
	}
	f += weights[from];

	return f;
}

graph* graph_from_file(const char* const filename)
{
	u8* const file = file_read_alloc(filename, alloc);

	char* line = (char*)file;
	char* eol = line;
	i8 size = 0;
	i8 capacity = 0;
	const char* name = NULL;
	i32 name_length = 0;
	graph* graph = NULL;

	// HEADER

	while (true)
	{
		while (*eol != '\n' && *eol != '\0') eol += 1;

		// graph name
		if (strncmp(NAME_PREFIX, line, sizeof(NAME_PREFIX) - 1) == 0)
		{
			while (*line != ':' && line != eol) line += 1;
			while (!isalnum(*line) && line != eol) line += 1;

			if (line != eol)
			{
				name = line;
				name_length = (i32)(eol - line);
			}
		}
		// graph size
		else if (strncmp(SIZE_PREFIX, line, sizeof(SIZE_PREFIX) - 1) == 0)
		{
			while (*line != ':' && line != eol) line += 1;

			if (line != eol)
			{
				line += 1;
				const i32 _size = strtol(line, NULL, 10);
				if (_size <= 0 || _size > 0x7F) abort();
				size = (i8)_size;
			}
		}
		// capacity
		else if (strncmp(CAPACITY_PREFIX, line, sizeof(CAPACITY_PREFIX) - 1) == 0)
		{
			while (*line != ':' && line != eol) line += 1;

			if (line != eol)
			{
				line += 1;
				const i32 _capacity = strtol(line, NULL, 10);
				if (_capacity <= 0 || _capacity > 0x7F) abort();
				capacity = (i8)_capacity;
			}
		}
		// header ending
		else if (strncmp(NODE_SECTION, line, sizeof(NODE_SECTION) - 1) == 0)
		{
			if (size == 0 || capacity == 0 || name == NULL) abort();
			graph = graph_init(size, name_length);
			graph->capacity = capacity;
			memcpy(graph->name, name, name_length);
			graph->name[name_length] = '\0';
			eol += 1;
			line = eol;
			break;
		}

		if (eol[0] == '\0' || eol[1] == '\0') abort();
		eol += 1;
		line = eol;
	}

	// NODES

	i8* xs = alloc_tmp(2 * size * sizeof(i8));
	i8* ys = xs + size;

	for (i8 i = 0; i < size; ++i)
	{
		while (*eol != '\n' && *eol != '\0') eol += 1;

		const i32 ind = strtol(line, &line, 10);
		const i32 x = strtol(line, &line, 10);
		const i32 y = strtol(line, &line, 10);

		if (ind < 1 || ind > size || x < 0 || x > 100 || y < 0 || y > 100) abort();

		xs[ind - 1] = (i8)x;
		ys[ind - 1] = (i8)y;

		if (eol[0] == '\0') abort();
		eol += 1;
		line = eol;
	}

	// CALCULATE WEIGHTS

	i16* res = graph->weights;
	for (i8 from = 0; from < size; ++from)
	{
		for (i8 to = 0; to < size; ++to)
		{
			const i8 dx = xs[to] - xs[from];
			const i8 dy = ys[to] - ys[from];

			const i32 d2 = dx * dx + dy * dy;
			const double d = sqrt(d2);

			const i16 scaled_dist = (i16)(d / SCALING * 0x7FFF);

			*res++ = scaled_dist;
		}
	}

	// DEMANDS

	if (strncmp(DEMAND_SECTION, line, sizeof(DEMAND_SECTION - 1)) != 0) abort();

	while (*eol != '\n' && *eol != '\0') eol += 1;
	if (eol[0] == '\0') abort();
	eol += 1;
	line = eol;

	for (i8 i = 0; i < size; ++i)
	{
		while (*eol != '\n' && *eol != '\0') eol += 1;

		const i32 ind = strtol(line, &line, 10);
		const i32 demand = strtol(line, &line, 10);

		if (ind < 1 || ind > size || demand < 0 || demand > capacity) abort();

		graph->demands[ind - 1] = (i8)demand;

		if (eol[0] == '\0') abort();
		eol += 1;
		line = eol;
	}

	dealloc(file);
	return graph;
}

void graph_deinit(graph* const graph)
{
	dealloc(graph);
}

i16 graph_get_weight(const graph* const graph, const i8 from, const i8 to)
{
	ASSERT(from >= 0 && to >= 0 && from < graph->size && to < graph->size);

	return graph->weights[from * graph->size + to];
}

i32 graph_solve_random(const graph* const graph, const i32 count)
{
	ASSERT(count >= 0);
	const i8 city_count = graph->size - 1;
	const i16 gene_count = graph->size * 2 - 3;

	i32 best_fitness = 0x7FFFFFFF;
	i8* const buf = alloc_tmp(gene_count * sizeof(i8));
	{
		i16 g;
		for (g = 0; g < city_count; ++g) buf[g] = (i8)(g + 1);
		for (; g < gene_count; ++g) buf[g] = 0;
	}

	for (i32 i = 0; i < count; ++i)
	{
		partial_shuffle_i8(buf, city_count - 1, city_count);

		i32 f = graph_fit_single(graph, buf);

		if (f < best_fitness)
		{
			best_fitness = f;
		}
	}

	return best_fitness;
}

i32 graph_solve_greedy(const graph* const graph, const i8 start_city)
{
	const i8* const demands = graph->demands;
	const i8 city_count = graph->size - 1;
	const i8 max_capacity = graph->capacity;

	i8* const buf = alloc_tmp(city_count * sizeof(i8));

	buf[0] = start_city;
	for (i8 i = 1; i < city_count; ++i) buf[i] = i + 1;
	buf[start_city - 1] = 1;

	i8 capacity_left = max_capacity - demands[start_city];
	ASSERT(capacity_left >= 0);
	i32 f = graph_get_weight(graph, 0, start_city);

	i8 from = start_city;

	for (i8 i = 0; i < city_count; ++i)
	{
		i8 shortest_to_ind = i;
		i16 shortest_length = graph_get_weight(graph, from, buf[shortest_to_ind]);
		for (i8 edge_ind = i + 1; edge_ind < city_count; ++edge_ind)
		{
			const i16 edge_length = graph_get_weight(graph, from, buf[edge_ind]);
			if (edge_length < shortest_length)
			{
				shortest_to_ind = edge_ind;
				shortest_length = edge_length;
			}
		}

		const i8 to = buf[shortest_to_ind];
		const i8 demand = demands[to];

		capacity_left -= demand;
		if (capacity_left < 0)
		{
			f += graph_get_weight(graph, from, 0);
			f += graph_get_weight(graph, 0, to);
			capacity_left = max_capacity - demand;
			ASSERT(capacity_left >= 0);
		}
		else
		{
			f += shortest_length;
		}

		const i8 tmp = buf[i];
		buf[i] = buf[shortest_to_ind];
		buf[shortest_to_ind] = tmp;

		from = to;
	}
	f += graph_get_weight(graph, from, 0);

	return f;
}

static graph* graph_init(const i8 size, const i32 name_length)
{
	ASSERT(size > 0 && name_length > 0);

	const i32 weights_offset = sizeof(graph);
	const i32 demands_offset = weights_offset + size * size * (i32)sizeof(i16);
	const i32 name_offset = demands_offset + size * sizeof(i8);
	const i32 data_size = name_offset + (name_length + 1) * (i32)sizeof(char);

	u8* mem = alloc(data_size);

	graph* const g = (graph*)mem;
	g->size = size;
	g->weights = (i16*)(mem + weights_offset);
	g->demands = (i8*)(mem + demands_offset);
	g->name = (char*)(mem + name_offset);

	return g;
}
