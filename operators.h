#pragma once

typedef struct graph graph;
typedef struct pop pop;
typedef struct tabu tabu;

// interface definitions

typedef void Iinit(const graph* graph, i8* genes, i32 count);

typedef struct Iselect
{
	i8* (*select)(const void* self, const pop* p);
	void (*preprocess)(const void* self, const pop* p);
	void (*select_data_init)(const void* self, pop* p);
	void (*select_data_deinit)(const void* self, pop* p);
} Iselect;

typedef void Icross(i16 gene_count, const i8* parent1, const i8* parent2, i8* child);

typedef void Imutate(i16 gene_count, i8* genes, i32 count, float chance);

typedef double Icooling(const void* self, double temp, i32 iter);

typedef void Ifithash(const tabu* t, i32* fit, u32* hash);

// Iinit objects

Iinit random_Iinit;
Iinit greedy_Iinit;

// Iselect objects

typedef struct tourney { float size; } tourney;
extern const Iselect tourney_Iselect;
extern const Iselect tourney_Iselect_distinct;

typedef struct roulette { float k; } roulette;
extern const Iselect roulette_Iselect;

// Icross objects

Icross ox_Icross;
Icross pmx_Icross;

// Imutation objects

Imutate swap_single_Imutate;
Imutate swap_many_Imutate;
Imutate invert_Imutate;
void mutate_copy(Imutate* op, i16 gene_count, const i8* src, i8* dst);

// Icooling objects

typedef struct cool { double k; double t_min; } cool;

Icooling cool_Icooling_subtract;
Icooling cool_Icooling_exp;

// Ifithash objects

Ifithash simple_Ifithash;
Ifithash complex_Ifithash;

void make_canonical_many(const graph* graph, i8* genes, i32 count);
void make_canonical_single(const graph* graph, i8* genes);
