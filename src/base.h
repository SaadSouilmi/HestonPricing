// Base layer: type aliases, keyword macros (internal/global/...), and
// assert helpers. Conventions taken in part from raddebugger's base layer
// (MIT, https://github.com/EpicGamesExt/raddebugger) and the Handmade Hero
// style (https://hero.handmade.network/).

#ifndef BASE_H
#define BASE_H

#include <stdint.h>

#define PI 3.14159265358979323846

////////////////////////////////////////////////////////////////
//~ Type Aliases

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef float f32;
typedef double f64;
typedef double _Complex f64c;

typedef i8 b8;
typedef i32 b32;

#define Min(a, b) (((a) < (b)) ? (a) : (b))
#define Max(a, b) (((a) > (b)) ? (a) : (b))
////////////////////////////////////////////////////////////////
//~ Build Configuration

#if !defined(BUILD_DEBUG)
#define BUILD_DEBUG 1
#endif

////////////////////////////////////////////////////////////////
//~ Keyword Macros

#define internal static
#define global static
#define local_persist static
#define read_only static const

////////////////////////////////////////////////////////////////
//~ Assert Macros

#if defined(_MSC_VER)
#define DEBUG_BREAK() __debugbreak()
#else
#define DEBUG_BREAK() __builtin_trap()
#endif

#define AssertAlways(expr)                                                     \
    do {                                                                       \
        if (!(expr)) {                                                         \
            DEBUG_BREAK();                                                     \
        }                                                                      \
    } while (0)

#if BUILD_DEBUG
#define Assert(expr) AssertAlways(expr)
#else
#define Assert(expr) ((void)0)
#endif

#define InvalidPath AssertAlways(!"Invalid path")
#define NotImplemented AssertAlways(!"Not implemented")

#endif
