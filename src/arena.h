// SPDX-License-Identifier: MIT
// Growable reserve/commit memory arena.
#ifndef ARENA_H
#define ARENA_H

#include "base.h"

#include <stdint.h>
#include <string.h>

////////////////////////////////////////////////////////////////
//~ Size helpers

#define KiB(n) ((u64)(n) << 10)
#define MiB(n) ((u64)(n) << 20)
#define GiB(n) ((u64)(n) << 30)

#define ALIGN_UP_POW2(x, p) (((u64)(x) + ((u64)(p) - 1)) & ~((u64)(p) - 1))

////////////////////////////////////////////////////////////////
//~ Arena
//
// A single virtual reservation that commits pages on demand as
// allocations grow. The arena's own header lives at the base of the
// reservation, so the first user allocation starts past it.

typedef struct mem_arena mem_arena;
struct mem_arena {
    u64 reserve_size; // total reserved bytes
    u64 commit_size;  // commit granularity
    u64 committed;    // bytes currently committed
    u64 pos;          // bytes used (offset from arena base)
};

#define ARENA_ALIGN (sizeof(void *))
#define ARENA_HEADER_SIZE ALIGN_UP_POW2(sizeof(mem_arena), ARENA_ALIGN)

internal mem_arena *arena_create(u64 reserve_size, u64 commit_size);
internal void arena_destroy(mem_arena *arena);
internal void *arena_push(mem_arena *arena, u64 size, b32 non_zero);
internal void arena_pop(mem_arena *arena, u64 size);
internal void arena_pop_to(mem_arena *arena, u64 pos);
internal void arena_clear(mem_arena *arena);

#define PUSH_STRUCT(arena, T) (T *)arena_push((arena), sizeof(T), 0)
#define PUSH_STRUCT_NZ(arena, T) (T *)arena_push((arena), sizeof(T), 1)
#define PUSH_ARRAY(arena, T, n)                                                \
    (T *)arena_push((arena), sizeof(T) * (u64)(n), 0)
#define PUSH_ARRAY_NZ(arena, T, n)                                             \
    (T *)arena_push((arena), sizeof(T) * (u64)(n), 1)

////////////////////////////////////////////////////////////////
//~ OS virtual memory layer

internal u32 os_page_size(void);
internal void *os_reserve(u64 size);
internal b32 os_commit(void *ptr, u64 size);
internal b32 os_release(void *ptr, u64 size);

#endif
