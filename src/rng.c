#include <math.h>
#include <string.h>

// xoshiro256** by David Blackman and Sebastiano Vigna (public domain, CC0)
//   https://prng.di.unimi.it/xoshiro256starstar.c
// State is seeded with splitmix64 (public domain), as recommended by the
// authors to spread a single 64-bit seed across the 256-bit state.

internal u64
rotl64(u64 x, i32 k) {
    return (x << k) | (x >> (64 - k));
}

internal u64
splitmix64(u64 *state) {
    u64 z = (*state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

internal void
ran_init(Ran *rng, u64 seed) {
    u64 sm = seed;
    rng->s[0] = splitmix64(&sm);
    rng->s[1] = splitmix64(&sm);
    rng->s[2] = splitmix64(&sm);
    rng->s[3] = splitmix64(&sm);
}

internal u64
ran_u64(Ran *rng) {
    u64 *s = rng->s;
    u64 result = rotl64(s[1] * 5, 7) * 9;
    u64 t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl64(s[3], 45);
    return result;
}

// Full-resolution uniform double in [0,1) via Allen Downey's method
//   "Generating Pseudo-random Floating-Point Values" (2007)
//   https://allendowney.com/research/rand/
internal f64
ran_f64_downey(Ran *rng) {
    u32 exponent = 0;
    u32 k = 0;
    do {
        k = clz_safe(ran_u64(rng));
        exponent += k;
    } while ((k == 64) && (exponent < 1022));
    exponent = 1022 - (exponent < 1022 ? exponent : 1022);

    if (!exponent) {
        return 0;
    }

    u64 mantissa = ran_u64(rng) & 0xFFFFFFFFFFFFF;
    if ((mantissa == 0) && (clz_safe(ran_u64(rng)) == 0)) {
        exponent++;
    }

    u64 bits = ((u64)exponent << 52) | mantissa;
    f64 result;
    memcpy(&result, &bits, sizeof(result));
    return result;
}

// Reference paper
// https://research.wu.ac.at/ws/files/18953249/document.pdf
internal u64
ran_poisson(Ran *rng, f64 lambda) {
    if (lambda <= 10) {
        f64 thresh = exp(-lambda);
        f64 p = 1.0;
        u64 k = 0;
        do {
            k++;
            p *= ran_f64_downey(rng);
        } while (p >= thresh);

        return k - 1;
    }

    f64 slambda = sqrt(lambda);
    f64 b = 0.931 + 2.53 * slambda;
    f64 a = -0.059 + 0.02483 * b;
    f64 inv_alpha = 1.1239 + 1.1328 / (b - 3.4);
    f64 v_r = 0.9277 - 3.6224 / (b - 2.0);

    do {
        f64 V = ran_f64_downey(rng);
        f64 U;
        if (V <= 0.86 * v_r) {
            U = V / v_r - 0.43;
            return (u64)(0.445 + lambda + U * (2 * a / (0.5 - fabs(U)) + b));
        }
        if (V >= v_r) {
            U = ran_f64_downey(rng) - 0.5;
        }
        if (V < v_r) {
            U = V / v_r - 0.93;
            U = Sign(U) * 0.5 - U;
            V = ran_f64_downey(rng) * v_r;
        }
        f64 us = 0.5 - fabs(U);
        if ((us < 0.013) && (V > us))
            continue;

        i64 k = (i64)(0.445 + lambda + U * (2.0 * a / us + b));
        V *= inv_alpha / (a / (us * us) + b);
        if ((k >= 10) &&
            (log(V * slambda) <=
             ((f64)k + 0.5) * log(lambda / (f64)k) - lambda -
                 0.5 * log(2.0 * PI) + (f64)k -
                 (1.0 / 12.0 - 1.0 / (360.0 * (f64)k * (f64)k)) / ((f64)k)))
            return k;
        if (k >= 0 && k <= 9 &&
            (log(V) <= k * log(lambda) - lambda - log_fact[k]))
            return k;
    } while (1);
}

// Code from numpy source
// https://github.com/numpy/numpy/blob/main/numpy/random/src/distributions/distributions.c
internal f64
ran_standard_exponential(Ran *rng) {
    u64 ri;
    u8 idx;
    f64 x;
    do {
        ri = ran_u64(rng);
        ri >>= 3;
        idx = ri & 0xFF;
        ri >>= 8;
        x = ri * we_f64[idx];
        if (ri < ke_f64[idx]) {
            return x; /* 98.9% of the time we return here 1st try */
        } else if (idx == 0) {
            return ziggurat_exp_r - log1p(-ran_f64_downey(rng));
        } else if ((fe_f64[idx - 1] - fe_f64[idx]) * ran_f64_downey(rng) +
                       fe_f64[idx] <
                   exp(-x)) {
            return x;
        }
    } while (1);
}

internal f64
ran_standard_normal(Ran *rng) {
    u64 r;
    i32 sign;
    u64 rabs;
    i32 idx;
    f64 x, xx, yy;
    for (;;) {
        /* r = e3n52sb8 */
        r = ran_u64(rng);
        idx = r & 0xff;
        r >>= 8;
        sign = r & 0x1;
        rabs = (r >> 1) & 0x000fffffffffffff;
        x = rabs * wi_f64[idx];
        if (sign & 0x1)
            x = -x;
        if (rabs < ki_f64[idx])
            return x; /* 99.3% of the time return here */
        if (idx == 0) {
            for (;;) {
                /* Switch to 1.0 - U to avoid log(0.0), see GH 13361 */
                xx = -ziggurat_nor_inv_r * log1p(-ran_f64_downey(rng));
                yy = -log1p(-ran_f64_downey(rng));
                if (yy + yy > xx * xx)
                    return ((rabs >> 8) & 0x1) ? -(ziggurat_nor_r + xx)
                                               : ziggurat_nor_r + xx;
            }
        } else {
            if (((fi_f64[idx - 1] - fi_f64[idx]) * ran_f64_downey(rng) +
                 fi_f64[idx]) < exp(-0.5 * x * x))
                return x;
        }
    }
}

internal f64
ran_gamma(Ran *rng, f64 shape) {
    f64 b, c;
    f64 U, V, X, Y;

    if (shape == 1.0) {
        return ran_standard_exponential(rng);
    } else if (shape == 0.0) {
        return 0.0;
    } else if (shape < 1.0) {
        for (;;) {
            U = ran_f64_downey(rng);
            V = ran_standard_exponential(rng);
            if (U <= 1.0 - shape) {
                X = pow(U, 1. / shape);
                if (X <= V) {
                    return X;
                }
            } else {
                Y = -log((1 - U) / shape);
                X = pow(1.0 - shape + shape * Y, 1. / shape);
                if (X <= (V + Y)) {
                    return X;
                }
            }
        }
    } else {
        b = shape - 1. / 3.;
        c = 1. / sqrt(9 * b);
        for (;;) {
            do {
                X = ran_standard_normal(rng);
                V = 1.0 + c * X;
            } while (V <= 0.0);

            V = V * V * V;
            U = ran_f64_downey(rng);
            if (U < 1.0 - 0.0331 * (X * X) * (X * X))
                return (b * V);
            /* log(0.0) ok here */
            if (log(U) < 0.5 * X * X + b * (1. - V + log(V)))
                return (b * V);
        }
    }
}

internal f64
ran_noncentral_chi2(Ran *rng, f64 d, f64 lambda) {
    if (d >= 1.0) {
        f64 z = ran_standard_normal(rng) + sqrt(lambda);
        return z * z + 2.0 * ran_gamma(rng, (d - 1.0) / 2.0);
    } else {
        u64 nn = ran_poisson(rng, lambda / 2.0);
        return 2.0 * ran_gamma(rng, (d + 2.0 * (f64)nn) / 2.0);
    }
}

////////////////////////////////////////////////////////////
//~ Helper functions

internal u32
clz_safe(u64 u) {
    return u ? __builtin_clzll(u) : 64;
}