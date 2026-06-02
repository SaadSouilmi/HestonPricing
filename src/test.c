#include <math.h>
#include <stdio.h>

#include "arena.h"
#include "heston.h"
#include "integrate.h"
#include "rng.h"

#include "arena.c"
#include "heston.c"
#include "integrate.c"
#include "rng.c"

////////////////////////////////////////////////////////////
//~ Reference Gauss-Legendre rule (N = 62) for validating
//~ compute_gauss_legendre_roots_and_weights.
//
// Source: standard tabulated Gauss-Legendre nodes/weights on [-1, 1].
// Reordered to ascending node order to match our solver's output
// convention (x[0] near -1, x[N-1] near +1). Weights are symmetric.

#define GL62_N 62

read_only f64 gl62_x[GL62_N] = {
    -0.9992598593087770, -0.9961022963162671, -0.9904299711892903,
    -0.9822559490972367, -0.9716007233716518, -0.9584911729739271,
    -0.9429604013923285, -0.9250476356362037, -0.9047981225210935,
    -0.8822630128318973, -0.8574992315120710, -0.8305693336040049,
    -0.8015413461039764, -0.7704885960554193, -0.7374895252831567,
    -0.7026274922222970, -0.6659905613354794, -0.6276712806468852,
    -0.5877664479530873, -0.5463768663002511, -0.5036070893447560,
    -0.4595651572401134, -0.4143623237171261, -0.3681127750465645,
    -0.3209333415941940, -0.2729432026967263, -0.2242635856041655,
    -0.1750174592490156, -0.1253292236158968, -0.0753243954962343,
    -0.0251292914218206, 0.0251292914218206,  0.0753243954962343,
    0.1253292236158968,  0.1750174592490156,  0.2242635856041655,
    0.2729432026967263,  0.3209333415941940,  0.3681127750465645,
    0.4143623237171261,  0.4595651572401134,  0.5036070893447560,
    0.5463768663002511,  0.5877664479530873,  0.6276712806468852,
    0.6659905613354794,  0.7026274922222970,  0.7374895252831567,
    0.7704885960554193,  0.8015413461039764,  0.8305693336040049,
    0.8574992315120710,  0.8822630128318973,  0.9047981225210935,
    0.9250476356362037,  0.9429604013923285,  0.9584911729739271,
    0.9716007233716518,  0.9822559490972367,  0.9904299711892903,
    0.9961022963162671,  0.9992598593087770,
};

read_only f64 gl62_w[GL62_N] = {
    0.0018992056795137, 0.0044163334569309, 0.0069260419018310,
    0.0094185794284204, 0.0118873901170105, 0.0143261918238065,
    0.0167288117901773, 0.0190891766585732, 0.0214013222776700,
    0.0236594072086828, 0.0258577269540247, 0.0279907281633146,
    0.0300530225739899, 0.0320394005816247, 0.0339448443794105,
    0.0357645406227681, 0.0374938925822800, 0.0391285317519631,
    0.0406643288824174, 0.0420974044103851, 0.0434241382580474,
    0.0446411789771244, 0.0457454522145702, 0.0467341684784155,
    0.0476048301841012, 0.0483552379634777, 0.0489834962205178,
    0.0494880179196993, 0.0498675285949524, 0.0501210695690433,
    0.0502480003752563, 0.0502480003752563, 0.0501210695690433,
    0.0498675285949524, 0.0494880179196993, 0.0489834962205178,
    0.0483552379634777, 0.0476048301841012, 0.0467341684784155,
    0.0457454522145702, 0.0446411789771244, 0.0434241382580474,
    0.0420974044103851, 0.0406643288824174, 0.0391285317519631,
    0.0374938925822800, 0.0357645406227681, 0.0339448443794105,
    0.0320394005816247, 0.0300530225739899, 0.0279907281633146,
    0.0258577269540247, 0.0236594072086828, 0.0214013222776700,
    0.0190891766585732, 0.0167288117901773, 0.0143261918238065,
    0.0118873901170105, 0.0094185794284204, 0.0069260419018310,
    0.0044163334569309, 0.0018992056795137,
};

////////////////////////////////////////////////////////////
//~ Test runner / logging

#define RunTest(call)                                                          \
    do {                                                                       \
        call;                                                                  \
        printf("  [ok] %s\n", #call);                                          \
    } while (0)

////////////////////////////////////////////////////////////
//~ Assertion macros for numerical tests

#define ExpectNear(actual, expected, tol)                                      \
    do {                                                                       \
        f64 _d = fabs((f64)((actual) - (expected)));                           \
        if (_d > (tol)) {                                                      \
            fprintf(stderr, "FAIL %s:%d  |%.17g - %.17g| = %.3g > %.3g\n",     \
                    __FILE__, __LINE__, (f64)(actual), (f64)(expected), _d,    \
                    (f64)(tol));                                               \
            DEBUG_BREAK();                                                     \
        }                                                                      \
    } while (0)

#define ExpectNearC(actual, expected, tol)                                     \
    do {                                                                       \
        f64c _a = (actual);                                                    \
        f64c _e = (expected);                                                  \
        f64 _d = cabs(_a - _e);                                                \
        if (_d > (tol)) {                                                      \
            fprintf(stderr,                                                    \
                    "FAIL %s:%d  |(%.17g+%.17gi) - (%.17g+%.17gi)| = %.3g\n",  \
                    __FILE__, __LINE__, creal(_a), cimag(_a), creal(_e),       \
                    cimag(_e), _d);                                            \
            DEBUG_BREAK();                                                     \
        }                                                                      \
    } while (0)

////////////////////////////////////////////////////////////
//~ Test cases — characteristic function

internal heston_params
default_params(void) {
    heston_params p = {
        .kappa = 1.5,
        .theta = 0.04,
        .sigma = 0.3,
        .rho = -0.7,
        .v0 = 0.04,
    };
    return p;
}

internal void
test_cf_normalisation(void) {
    /* phi_T(0) == 1 */
    heston_params p = default_params();
    f64c v = heston_cf(&p, 0.0 + 0.0 * I, 1.0);
    ExpectNearC(v, 1.0 + 0.0 * I, 1e-12);
}

internal void
test_cf_martingale(void) {
    /* phi_T(-i) == 1  (martingale property of e^{X_T}) */
    heston_params p = default_params();
    f64c v = heston_cf(&p, 0.0 - 1.0 * I, 1.0);
    ExpectNearC(v, 1.0 + 0.0 * I, 1e-12);
}

internal void
test_cf_conjugate_symmetry(void) {
    /* For real u: phi_T(-u) == conj(phi_T(u)) */
    heston_params p = default_params();
    f64c a = heston_cf(&p, 1.3 + 0.0 * I, 1.0);
    f64c b = heston_cf(&p, -1.3 + 0.0 * I, 1.0);
    ExpectNearC(b, conj(a), 1e-12);
}

internal void
test_gl_roots_and_weights_62(mem_arena *arena) {
    gl_quad q = {.N = 62};
    q.x = PUSH_ARRAY(arena, f64, q.N);
    q.w = PUSH_ARRAY(arena, f64, q.N);
    compute_gauss_legendre_roots_and_weights(&q);
    for (u16 i = 0; i < q.N; i++) {
        ExpectNear(q.x[i], gl62_x[i], 1e-13);
        ExpectNear(q.w[i], gl62_w[i], 1e-13);
    }
    arena_clear(arena);
}

internal void
test_poisson_moments(Ran *rng, f64 lambda) {
    u64 n = 2000000;
    f64 sum = 0.0, sum2 = 0.0;
    u64 zeros = 0;
    for (u64 i = 0; i < n; i++) {
        u64 k = ran_poisson(rng, lambda);
        sum += (f64)k;
        sum2 += (f64)k * (f64)k;
        if (k == 0)
            zeros++;
    }
    f64 mean = sum / (f64)n;
    f64 var = sum2 / (f64)n - mean * mean;
    f64 p0 = (f64)zeros / (f64)n;

    /* 5-sigma bands. SE(mean) = sqrt(lambda/n);
     * SE(var) = sqrt((mu4 - mu2^2)/n) with mu2=lambda, mu4=lambda+3lambda^2. */
    f64 mean_tol = 5.0 * sqrt(lambda / (f64)n);
    f64 var_tol = 5.0 * sqrt((lambda + 2.0 * lambda * lambda) / (f64)n);
    f64 p0_se = sqrt(exp(-lambda) * (1.0 - exp(-lambda)) / (f64)n);

    ExpectNear(mean, lambda, mean_tol);               /* E[k]   = lambda */
    ExpectNear(var, lambda, var_tol);                 /* Var[k] = lambda */
    ExpectNear(p0, exp(-lambda), 5.0 * p0_se + 1e-9); /* P(0) = e^-lambda */
}

internal void
test_poisson_chi2(mem_arena *arena, Ran *rng, f64 lambda) {
    /* Chi-squared goodness-of-fit against the exact Poisson PMF.
     * Catches shape errors that the moment test misses. */
    u64 n = 2000000;
    u32 kmax = (u32)(lambda + 12.0 * sqrt(lambda)) + 20; /* generous support */
    u64 *obs = PUSH_ARRAY(arena, u64, kmax + 1);         /* zeroed by arena */

    for (u64 i = 0; i < n; i++) {
        u64 k = ran_poisson(rng, lambda);
        if (k > kmax)
            k = kmax; /* lump far tail */
        obs[k]++;
    }

    /* Adaptive binning: accumulate until expected >= 5, then close the bin.
     * pmf iterates as p(k+1) = p(k) * lambda / (k+1), p(0) = e^-lambda. */
    f64 pmf = exp(-lambda);
    f64 chi2 = 0.0, bin_o = 0.0, bin_e = 0.0;
    u32 nbins = 0;
    for (u32 k = 0; k <= kmax; k++) {
        bin_o += (f64)obs[k];
        bin_e += (f64)n * pmf;
        if (bin_e >= 5.0 && k < kmax) {
            f64 d = bin_o - bin_e;
            chi2 += d * d / bin_e;
            nbins++;
            bin_o = 0.0;
            bin_e = 0.0;
        }
        pmf *= lambda / (f64)(k + 1);
    }
    if (bin_e > 0.0) { /* final bin + tail */
        f64 d = bin_o - bin_e;
        chi2 += d * d / bin_e;
        nbins++;
    }

    /* chi2 has mean = dof, std = sqrt(2 dof). 5-sigma band => no false fails,
     * but a shape bug pushes chi2 far past it. */
    u32 dof = nbins - 1;
    f64 band = 5.0 * sqrt(2.0 * (f64)dof);
    ExpectNear(chi2, (f64)dof, band);

    arena_clear(arena);
}

internal f64
std_normal_cdf(f64 x) {
    return 0.5 * erfc(-x * 0.70710678118654752440); /* Phi(x), 1/sqrt(2) */
}

internal void
test_normal_moments(Ran *rng) {
    /* Standard normal raw moments: E[X]=0, E[X^2]=1, E[X^3]=0, E[X^4]=3.
     * The 4th moment (kurtosis) is the real shape check. */
    u64 n = 5000000;
    f64 s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0;
    for (u64 i = 0; i < n; i++) {
        f64 x = ran_standard_normal(rng);
        f64 x2 = x * x;
        s1 += x;
        s2 += x2;
        s3 += x2 * x;
        s4 += x2 * x2;
    }
    f64 m1 = s1 / (f64)n, m2 = s2 / (f64)n, m3 = s3 / (f64)n, m4 = s4 / (f64)n;

    ExpectNear(m1, 0.0, 0.01); /* mean      */
    ExpectNear(m2, 1.0, 0.01); /* variance  */
    ExpectNear(m3, 0.0, 0.02); /* skew num. */
    ExpectNear(m4, 3.0, 0.05); /* kurtosis  */
}

internal void
test_normal_chi2(mem_arena *arena, Ran *rng) {
    /* Chi-squared goodness-of-fit against Phi. Bins over [-8, 8]; samples
     * outside clamp into the edge cells, which carry the tail mass. */
    u64 n = 2000000;
    f64 lo = -8.0, step = 0.05;
    u32 ncells = 320; /* (8 - -8) / 0.05 */
    u64 *obs = PUSH_ARRAY(arena, u64, ncells);

    for (u64 i = 0; i < n; i++) {
        f64 x = ran_standard_normal(rng);
        i32 idx = (i32)((x - lo) / step);
        if (idx < 0)
            idx = 0;
        if (idx >= (i32)ncells)
            idx = (i32)ncells - 1;
        obs[idx]++;
    }

    /* Adaptive binning: cell 0 = (-inf, lo+step), cell last = [hi-step, +inf).
     */
    f64 chi2 = 0.0, bin_o = 0.0, bin_e = 0.0, prev_cdf = 0.0;
    u32 nbins = 0;
    for (u32 j = 0; j < ncells; j++) {
        f64 cdf =
            (j == ncells - 1) ? 1.0 : std_normal_cdf(lo + (f64)(j + 1) * step);
        f64 cell_e = (f64)n * (cdf - prev_cdf);
        prev_cdf = cdf;

        bin_o += (f64)obs[j];
        bin_e += cell_e;
        if (bin_e >= 5.0 && j < ncells - 1) {
            f64 d = bin_o - bin_e;
            chi2 += d * d / bin_e;
            nbins++;
            bin_o = 0.0;
            bin_e = 0.0;
        }
    }
    if (bin_e > 0.0) {
        f64 d = bin_o - bin_e;
        chi2 += d * d / bin_e;
        nbins++;
    }

    u32 dof = nbins - 1;
    f64 band = 5.0 * sqrt(2.0 * (f64)dof);
    ExpectNear(chi2, (f64)dof, band);

    arena_clear(arena);
}

internal void
test_exponential_moments(Ran *rng) {
    /* Exp(1): mean = 1, var = 1; survival P(X>x) = e^-x. */
    u64 n = 2000000;
    f64 sum = 0.0, sum2 = 0.0;
    u64 gt1 = 0, gt2 = 0, gt3 = 0;
    for (u64 i = 0; i < n; i++) {
        f64 x = ran_standard_exponential(rng);
        sum += x;
        sum2 += x * x;
        if (x > 1.0)
            gt1++;
        if (x > 2.0)
            gt2++;
        if (x > 3.0)
            gt3++;
    }
    f64 mean = sum / (f64)n;
    f64 var = sum2 / (f64)n - mean * mean;

    ExpectNear(mean, 1.0, 0.01);                     /* E[X]   = 1 */
    ExpectNear(var, 1.0, 0.02);                      /* Var[X] = 1 */
    ExpectNear((f64)gt1 / (f64)n, exp(-1.0), 0.005); /* P(X>1) */
    ExpectNear((f64)gt2 / (f64)n, exp(-2.0), 0.005); /* P(X>2) */
    ExpectNear((f64)gt3 / (f64)n, exp(-3.0), 0.005); /* P(X>3) */
}

internal void
test_exponential_chi2(mem_arena *arena, Ran *rng) {
    /* Chi-squared goodness-of-fit against Exp(1). Bins over [0, 20];
     * cell expected = n*(e^-a - e^-b) via the survival function. */
    u64 n = 2000000;
    f64 step = 0.05;
    u32 ncells = 400; /* 20 / 0.05 */
    u64 *obs = PUSH_ARRAY(arena, u64, ncells);

    for (u64 i = 0; i < n; i++) {
        f64 x = ran_standard_exponential(rng);
        i32 idx = (i32)(x / step);
        if (idx < 0)
            idx = 0;
        if (idx >= (i32)ncells)
            idx = (i32)ncells - 1;
        obs[idx]++;
    }

    f64 chi2 = 0.0, bin_o = 0.0, bin_e = 0.0, prev_surv = 1.0; /* S(0) = 1 */
    u32 nbins = 0;
    for (u32 j = 0; j < ncells; j++) {
        f64 surv = (j == ncells - 1) ? 0.0 : exp(-(f64)(j + 1) * step);
        f64 cell_e = (f64)n * (prev_surv - surv);
        prev_surv = surv;

        bin_o += (f64)obs[j];
        bin_e += cell_e;
        if (bin_e >= 5.0 && j < ncells - 1) {
            f64 d = bin_o - bin_e;
            chi2 += d * d / bin_e;
            nbins++;
            bin_o = 0.0;
            bin_e = 0.0;
        }
    }
    if (bin_e > 0.0) {
        f64 d = bin_o - bin_e;
        chi2 += d * d / bin_e;
        nbins++;
    }

    u32 dof = nbins - 1;
    f64 band = 5.0 * sqrt(2.0 * (f64)dof);
    ExpectNear(chi2, (f64)dof, band);

    arena_clear(arena);
}

internal f64
gamma_p(f64 a, f64 x) {
    /* Regularized lower incomplete gamma P(a,x), a>0, x>=0.
     * Standard evaluation from DLMF (https://dlmf.nist.gov/8):
     *   - power series for the lower function, 8.11.4, when x < a+1;
     *   - continued fraction for the upper function Q, 8.9.2, otherwise,
     *     summed with the modified Lentz algorithm (Lentz 1976).
     * Cf. the Cephes library (BSD-licensed) igam/igamc for a reference impl. */
    if (x <= 0.0)
        return 0.0;
    f64 log_norm = a * log(x) - x - lgamma(a); /* log(x^a e^-x / Gamma(a)) */

    if (x < a + 1.0) {
        /* DLMF 8.11.4:  P = log_norm-factor * sum_{k>=0} x^k / (a)_{k+1}. */
        f64 term = 1.0 / a;
        f64 sum = term;
        f64 denom = a;
        for (i32 k = 0; k < 1000; k++) {
            denom += 1.0;
            term *= x / denom;
            sum += term;
            if (fabs(term) < fabs(sum) * 1e-16)
                break;
        }
        return sum * exp(log_norm);
    } else {
        /* DLMF 8.9.2 continued fraction for Q(a,x), modified Lentz. */
        f64 tiny = 1e-300;
        f64 bj = x + 1.0 - a;
        f64 cj = 1.0 / tiny;
        f64 dj = 1.0 / bj;
        f64 frac = dj;
        for (i32 j = 1; j < 1000; j++) {
            f64 aj = -(f64)j * ((f64)j - a);
            bj += 2.0;
            dj = aj * dj + bj;
            if (fabs(dj) < tiny)
                dj = tiny;
            cj = bj + aj / cj;
            if (fabs(cj) < tiny)
                cj = tiny;
            dj = 1.0 / dj;
            f64 ratio = dj * cj;
            frac *= ratio;
            if (fabs(ratio - 1.0) < 1e-16)
                break;
        }
        f64 q = exp(log_norm) * frac;
        return 1.0 - q;
    }
}

internal void
test_gamma_moments(Ran *rng, f64 shape) {
    /* Gamma(shape, 1): mean = shape, var = shape. */
    u64 n = 2000000;
    f64 sum = 0.0, sum2 = 0.0;
    for (u64 i = 0; i < n; i++) {
        f64 x = ran_gamma(rng, shape);
        sum += x;
        sum2 += x * x;
    }
    f64 mean = sum / (f64)n;
    f64 var = sum2 / (f64)n - mean * mean;

    ExpectNear(mean, shape, 0.02 * shape + 0.005);
    ExpectNear(var, shape, 0.05 * shape + 0.01);
}

internal void
test_gamma_chi2(mem_arena *arena, Ran *rng, f64 shape) {
    /* Chi-squared goodness-of-fit against Gamma(shape, 1) via gamma_p. */
    u64 n = 2000000;
    f64 xmax = shape + 15.0 * sqrt(shape) + 20.0;
    u32 ncells = 500;
    f64 step = xmax / (f64)ncells;
    u64 *obs = PUSH_ARRAY(arena, u64, ncells);

    for (u64 i = 0; i < n; i++) {
        f64 x = ran_gamma(rng, shape);
        i32 idx = (i32)(x / step);
        if (idx < 0)
            idx = 0;
        if (idx >= (i32)ncells)
            idx = (i32)ncells - 1;
        obs[idx]++;
    }

    f64 chi2 = 0.0, bin_o = 0.0, bin_e = 0.0, prev_cdf = 0.0; /* P(a,0)=0 */
    u32 nbins = 0;
    for (u32 j = 0; j < ncells; j++) {
        f64 cdf = (j == ncells - 1) ? 1.0 : gamma_p(shape, (f64)(j + 1) * step);
        f64 cell_e = (f64)n * (cdf - prev_cdf);
        prev_cdf = cdf;

        bin_o += (f64)obs[j];
        bin_e += cell_e;
        if (bin_e >= 5.0 && j < ncells - 1) {
            f64 d = bin_o - bin_e;
            chi2 += d * d / bin_e;
            nbins++;
            bin_o = 0.0;
            bin_e = 0.0;
        }
    }
    if (bin_e > 0.0) {
        f64 d = bin_o - bin_e;
        chi2 += d * d / bin_e;
        nbins++;
    }

    u32 dof = nbins - 1;
    f64 band = 5.0 * sqrt(2.0 * (f64)dof);
    ExpectNear(chi2, (f64)dof, band);

    arena_clear(arena);
}

int
main(void) {
    mem_arena *arena = arena_create(GiB(1), MiB(1));

    RunTest(test_cf_normalisation());
    RunTest(test_cf_martingale());
    RunTest(test_cf_conjugate_symmetry());
    RunTest(test_gl_roots_and_weights_62(arena));

    Ran rng;
    ran_init(&rng, 0x9e3779b97f4a7c15ULL);
    RunTest(test_poisson_moments(&rng, 3.0));      /* Knuth branch */
    RunTest(test_poisson_moments(&rng, 50.0));     /* PTRD branch  */
    RunTest(test_poisson_chi2(arena, &rng, 3.0));  /* Knuth branch */
    RunTest(test_poisson_chi2(arena, &rng, 50.0)); /* PTRD branch  */
    RunTest(test_normal_moments(&rng));
    RunTest(test_normal_chi2(arena, &rng));
    RunTest(test_exponential_moments(&rng));
    RunTest(test_exponential_chi2(arena, &rng));
    RunTest(test_gamma_moments(&rng, 0.5));     /* GS branch (shape<1) */
    RunTest(test_gamma_moments(&rng, 2.5));     /* MT branch (shape>=1) */
    RunTest(test_gamma_chi2(arena, &rng, 0.5)); /* GS branch */
    RunTest(test_gamma_chi2(arena, &rng, 2.5)); /* MT branch */

    printf("all tests passed\n");

    arena_destroy(arena);
    return 0;
}
