#include "stdafx.h"
#include "tester.h"
#include "graph.h"
#include "metah.h"

void test_metah(FILE* const file, const void* const config, const Imetah* const metah_impl, const i32 iterations, const i32 repeats)
{
	ASSERT(iterations > 0 && repeats > 0);

	i32 best = 0x7FFFFFFF;
	i32 worst = 0;

	i64 sum = 0;
	i64 sum2 = 0;

	void* const metah = metah_impl->init(config);

	for (i32 repeat = 0; repeat < repeats; ++repeat)
	{
		metah_impl->reset(metah);

		i32 f = -1;
		for (i32 i = 0; i < iterations; ++i)
		{
			f = metah_impl->iterate(metah);
		}

		printf("\33[2K\r");
		//metah_impl->fprint_overall_best(metah, stdout);
		printf("[%d] %d", repeat + 1, convert_units(f));

		if (f < best) best = f;
		if (f > worst) worst = f;
		sum += f;
		sum2 += (i64)f * f;
	}
	printf("\33[2K\r");

	metah_impl->deinit(metah);

	const double inv_repeats = 1.0 / repeats;
	const double avg = sum * inv_repeats;
	const double inv_repeats_m1 = 1.0 / (repeats - 1.0);
	const double std = repeats >= 2 ? sqrt(sum2 * inv_repeats_m1 - sum * avg * inv_repeats_m1) : 0.0;

	fprintf(file, "%d\t%d\t%d\t%d\n", convert_units(best), convert_units(worst), convert_units(avg), convert_units(std));
}

void test_random(FILE* const file, const graph* const graph, const i32 count, const i32 repeats)
{
	ASSERT(count > 0 && repeats > 0);

	i32 best = 0x7FFFFFFF;
	i32 worst = 0;

	i64 sum = 0;
	i64 sum2 = 0;

	for (i32 repeat = 0; repeat < repeats; ++repeat)
	{
		const i32 f = graph_solve_random(graph, count);
		printf("\r[%d] %d", repeat + 1, convert_units(f));

		if (f < best) best = f;
		if (f > worst) worst = f;
		sum += f;
		sum2 += (i64)f * f;
	}
	printf("\33[2K\r");

	const double inv_repeats = 1.0 / repeats;
	const double avg = sum * inv_repeats;
	const double inv_repeats_m1 = 1.0 / (repeats - 1.0);
	const double std = repeats >= 2 ? sqrt(sum2 * inv_repeats_m1 - sum * avg * inv_repeats_m1) : 0.0;

	fprintf(file, "%d\t%d\t%d\t%d\n", convert_units(best), convert_units(worst), convert_units(avg), convert_units(std));
}

void test_greedy(FILE* const file, const graph* const graph)
{
	const i8 city_count = graph->size - 1;

	i32 best = 0x7FFFFFFF;
	i32 worst = 0;

	i64 sum = 0;
	i64 sum2 = 0;

	for (i8 start_city = 1; start_city <= city_count; ++start_city)
	{
		const i32 f = graph_solve_greedy(graph, start_city);
		printf("\r[%d] %d", start_city, convert_units(f));

		if (f < best) best = f;
		if (f > worst) worst = f;
		sum += f;
		sum2 += (i64)f * f;
	}
	printf("\33[2K\r");

	const double inv_repeats = 1.0 / city_count;
	const double avg = sum * inv_repeats;
	const double inv_repeats_m1 = 1.0 / (city_count - 1.0);
	const double std = city_count >= 2 ? sqrt(sum2 * inv_repeats_m1 - sum * avg * inv_repeats_m1) : 0.0;

	fprintf(file, "%d\t%d\t%d\t%d\n", convert_units(best), convert_units(worst), convert_units(avg), convert_units(std));
}

void test_detail(const void* const config, const Imetah* const metah_impl, const i32 iterations, const char* const filename)
{
	ASSERT(iterations > 0);

	FILE* file;
	fopen_s(&file, filename, "wb");
	if (!file) abort();

	void* const metah = metah_impl->init(config);
	metah_impl->reset(metah);
	metah_impl->fprint_detail_header(metah, file);
	metah_impl->fprint_detail(metah, file);

	for (i32 i = 0; i < iterations; ++i)
	{
		metah_impl->iterate(metah);
		metah_impl->fprint_detail(metah, file);
	}

	metah_impl->deinit(metah);
	fclose(file);
}
