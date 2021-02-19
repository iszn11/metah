/* Based on implementation from https://www.pcg-random.org
 *
 * Original licensed under the Apache License, Version 2.0 (the "License")
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * A copy of the License is available at http://www.apache.org/licenses/LICENSE-2.0
 */

#include "stdafx.h"

static rng_state rng;

void random_seed(u64 seed, u64 seq)
{
	rng.state = 0;
	rng.inc = (seq << 1U) | 1U;
	random_u32();
	rng.state += seed;
}

void random_set_state(const rng_state* const state)
{
	rng.state = state->state;
	rng.inc = state->inc;
}

rng_state* random_get_state()
{
	return &rng;
}

u32 random_u32()
{
	const u64 old_state = rng.state;
	rng.state = old_state * 6364136223846793005ULL + rng.inc;
	const u32 xorshifted = (u32)(((old_state >> 18U) ^ old_state) >> 27U);
	const u32 rot = old_state >> 59U;
	#pragma warning(suppress : 4146) // intended unsigned negation
	return (xorshifted >> rot) | (xorshifted << ((-rot) & 0x1FU));
}

i32 random_range(i32 a, i32 b)
{
	ASSERT(a < b);
	const u32 diff = (u32)b - (u32)a;
	#pragma warning(suppress : 4146) // intended unsigned negation
	const u32 threshold = -diff % diff;

	while (true)
	{
		const u32 r = random_u32();
		if (r >= threshold)
		{
			return r % diff + (u32)a;
		}
	}
}

float random_float()
{
	return (float)ldexp(random_u32(), -32);
}

bool random_chance(float chance)
{
	const double r = ldexp(random_u32(), -32);
	return r < chance;
}

void partial_shuffle_i8(i8* const array, const i32 count, const i32 size)
{
	ASSERT(count > 0 && count < size);

	for (i32 i = 0; i < count; ++i)
	{
		const i32 rand_ind = random_range(i, size);

		const i8 tmp = array[i];
		array[i] = array[rand_ind];
		array[rand_ind] = tmp;
	}
}

void partial_shuffle_i16(i16* const array, const i32 count, const i32 size)
{
	ASSERT(count > 0 && count < size);

	for (i16 i = 0; i < count; ++i)
	{
		const i32 rand_ind = random_range(i, size);

		const i16 tmp = array[i];
		array[i] = array[rand_ind];
		array[rand_ind] = tmp;
	}
}
