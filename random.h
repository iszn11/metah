/* Based on implementation from https://www.pcg-random.org
 *
 * Original licensed under the Apache License, Version 2.0 (the "License")
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * A copy of the License is available at http://www.apache.org/licenses/LICENSE-2.0
 */

#pragma once

typedef struct rng_state
{
	u64 state;
	u64 inc;
} rng_state;

void random_seed(u64 seed, u64 seq);
rng_state* random_get_state();

u32 random_u32();
i32 random_range(i32 a, i32 b);
float random_float();
bool random_chance(float chance);

void partial_shuffle_i8(i8* array, i32 count, i32 size);
void partial_shuffle_i16(i16* array, i32 count, i32 size);
