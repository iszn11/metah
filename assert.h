#pragma once

#ifndef NDEBUG
#define ASSERT(cond) do { if (!(cond)) __debugbreak(); } while (0)
#else
#define ASSERT(cond)
#endif
