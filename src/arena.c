// SPDX-License-Identifier: MIT
//
// Growable reserve/commit arena. One large virtual reservation is made up
// front; physical pages are committed lazily in `commit_size` chunks as the
// bump pointer advances. The mem_arena header is stored at the base of the
// reservation itself, so no separate allocation is needed for it.

#if defined(__linux__) || defined(__APPLE__)
#define _DEFAULT_SOURCE
#include <sys/mman.h>
#include <unistd.h>
#elif defined(_WIN32)
#include <windows.h>
#endif

#include "arena.h"

internal mem_arena *
arena_create(u64 reserve_size, u64 commit_size) {
    u32 page = os_page_size();
    reserve_size = ALIGN_UP_POW2(reserve_size, page);
    commit_size = ALIGN_UP_POW2(commit_size, page);

    u8 *base = (u8 *)os_reserve(reserve_size);
    if (base == 0) {
        return 0;
    }

    u64 first_commit = ALIGN_UP_POW2(ARENA_HEADER_SIZE, commit_size);
    if (!os_commit(base, first_commit)) {
        os_release(base, reserve_size);
        return 0;
    }

    mem_arena *arena = (mem_arena *)base;
    arena->reserve_size = reserve_size;
    arena->commit_size = commit_size;
    arena->committed = first_commit;
    arena->pos = ARENA_HEADER_SIZE;
    return arena;
}

internal void
arena_destroy(mem_arena *arena) {
    os_release(arena, arena->reserve_size);
}

internal void *
arena_push(mem_arena *arena, u64 size, b32 non_zero) {
    u64 start = ALIGN_UP_POW2(arena->pos, ARENA_ALIGN);
    u64 end = start + size;
    if (end > arena->reserve_size) {
        return 0; // exhausted the reservation
    }

    if (end > arena->committed) {
        u64 target = ALIGN_UP_POW2(end, arena->commit_size);
        if (target > arena->reserve_size) {
            target = arena->reserve_size;
        }
        u8 *from = (u8 *)arena + arena->committed;
        if (!os_commit(from, target - arena->committed)) {
            return 0;
        }
        arena->committed = target;
    }

    arena->pos = end;
    u8 *out = (u8 *)arena + start;
    if (!non_zero) {
        memset(out, 0, size);
    }
    return out;
}

internal void
arena_pop(mem_arena *arena, u64 size) {
    u64 floor = ARENA_HEADER_SIZE;
    u64 available = arena->pos - floor;
    if (size > available) {
        size = available;
    }
    arena->pos -= size;
}

internal void
arena_pop_to(mem_arena *arena, u64 pos) {
    if (pos < arena->pos) {
        arena_pop(arena, arena->pos - pos);
    }
}

internal void
arena_clear(mem_arena *arena) {
    arena_pop_to(arena, ARENA_HEADER_SIZE);
}

////////////////////////////////////////////////////////////////
//~ OS virtual memory layer

#if defined(__linux__) || defined(__APPLE__)

internal u32
os_page_size(void) {
    return (u32)sysconf(_SC_PAGESIZE);
}

internal void *
os_reserve(u64 size) {
    void *p = mmap(0, size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    return (p == MAP_FAILED) ? 0 : p;
}

internal b32
os_commit(void *ptr, u64 size) {
    return mprotect(ptr, size, PROT_READ | PROT_WRITE) == 0;
}

internal b32
os_release(void *ptr, u64 size) {
    return munmap(ptr, size) == 0;
}

#elif defined(_WIN32)

internal u32
os_page_size(void) {
    SYSTEM_INFO info = {0};
    GetSystemInfo(&info);
    return info.dwPageSize;
}

internal void *
os_reserve(u64 size) {
    return VirtualAlloc(0, size, MEM_RESERVE, PAGE_READWRITE);
}

internal b32
os_commit(void *ptr, u64 size) {
    return VirtualAlloc(ptr, size, MEM_COMMIT, PAGE_READWRITE) != 0;
}

internal b32
os_release(void *ptr, u64 size) {
    (void)size; // VirtualFree(MEM_RELEASE) requires size 0
    return VirtualFree(ptr, 0, MEM_RELEASE);
}

#endif
