#pragma once

#define GRAPHS_COUNT 7

typedef struct graph graph;

extern graph* GRAPHS[GRAPHS_COUNT];

void graphs_init();
void graphs_deinit();
