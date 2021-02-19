#include "stdafx.h"

#include "graph.h"
#include "graphs.h"

#include "tester.h"
#include "pop.h"
#include "tabu.h"
#include "sa.h"

#include <time.h>

const i32 REPEATS = 1;

i32 main()
{
	// SEED

	const u64 seed = time(NULL);
	const u64 seq = 0;

	printf("SEED: %lld,%lld\n", seed, seq);
	random_seed(seed, seq);

	// GRAPH INITIALIZATION

	graphs_init();

	// WORK

	const tourney tourney = { .size = 0.01f };
	const void* const tourney_ptr = &tourney;
	pop_config p_cfg = {
		.init_fn = greedy_Iinit,
		.select = tourney_ptr,
		.select_impl = &tourney_Iselect,
		.cross_fn = pmx_Icross,
		.mutate_fn = invert_Imutate,

		.cross_chance = 0.0f,
		.mutate_chance = 0.5f,
		.force_canonical = false,

		.pop_size = 500
	};
	const i32 p_iter = 500;

	tabu_config t_cfg = {
		.init_fn = greedy_Iinit,
		.neighbor_fn = invert_Imutate,
		.fithash_fn = complex_Ifithash,

		.tabu_capacity = 1,
		.neighbor_count = 200
	};
	const i32 t_iter = 4000;

	const cool cool = { .k = 0.9996, .t_min = 1000.0 };
	const void* const cool_ptr = &cool;
	sa_config s_cfg = {
		.init_fn = greedy_Iinit,
		.neighbor_fn = invert_Imutate,
		.cooling = cool_ptr,
		.cooling_impl = cool_Icooling_exp,

		.initial_temp = 50000.0,
		.neighbor_count = 4
	};
	const i32 s_iter = 200000;

	const i32 random_count = 800000;

	for (i32 i = 0; i < GRAPHS_COUNT; ++i)
	{
		printf("=== [%d] %s ===\n", i + 1, GRAPHS[i]->name);

		p_cfg.graph = GRAPHS[i];
		t_cfg.graph = GRAPHS[i];
		s_cfg.graph = GRAPHS[i];

		test_metah(stdout, &p_cfg, &pop_Imetah, p_iter, REPEATS);
		test_metah(stdout, &t_cfg, &tabu_Imetah, t_iter, REPEATS);
		test_metah(stdout, &s_cfg, &sa_Imetah, s_iter, REPEATS);

		test_random(stdout, GRAPHS[i], random_count, REPEATS);
		test_greedy(stdout, GRAPHS[i]);
	}

	// DEINITIALIZATION

	graphs_deinit();
}
