#pragma once

typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned int u32;
typedef unsigned long long u64;

typedef signed char i8;
typedef short i16;
typedef int i32;
typedef long long i64;

_STATIC_ASSERT(sizeof(u8) == 1);
_STATIC_ASSERT(sizeof(u16) == 2);
_STATIC_ASSERT(sizeof(u32) == 4);
_STATIC_ASSERT(sizeof(u64) == 8);

_STATIC_ASSERT(sizeof(i8) == 1);
_STATIC_ASSERT(sizeof(i16) == 2);
_STATIC_ASSERT(sizeof(i32) == 4);
_STATIC_ASSERT(sizeof(i64) == 8);
