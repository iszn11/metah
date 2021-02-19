#pragma once

typedef struct graph
{
	i8 size;
	i8 capacity;
	i16* weights;
	i8* demands;
	char* name;
} graph;

i32 convert_units(double internal_dist);

void graph_fit_many(const graph* graph, const i8* genes, i32* fitness, i32 count, i32* best_ind, i32* best_fit, i32* worst_fit);
i32 graph_fit_single(const graph* graph, const i8* genes);

graph* graph_from_file(const char* filename);
void graph_deinit(graph* graph);

i16 graph_get_weight(const graph* graph, i8 from, i8 to);

i32 graph_solve_random(const graph* graph, i32 count);
i32 graph_solve_greedy(const graph* graph, i8 start_city);
