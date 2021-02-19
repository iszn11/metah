#include "stdafx.h"
#include "graphs.h"
#include "graph.h"

static const char* GRAPH_FILENAMES[GRAPHS_COUNT] =
{
	"A-n32-k5.vrp",
	"A-n37-k6.vrp",
	"A-n39-k5.vrp",
	"A-n45-k6.vrp",
	"A-n48-k7.vrp",
	"A-n54-k7.vrp",
	"A-n60-k9.vrp",
};

graph* GRAPHS[GRAPHS_COUNT];

void graphs_init()
{
	for (i32 i = 0; i < GRAPHS_COUNT; ++i)
	{
		GRAPHS[i] = graph_from_file(GRAPH_FILENAMES[i]);
	}
}

void graphs_deinit()
{
	for (i32 i = 0; i < GRAPHS_COUNT; ++i)
	{
		graph_deinit(GRAPHS[i]);
	}
}
