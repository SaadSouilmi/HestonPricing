#ifndef HESTON_H
#define HESTON_H

#include "base.h"
#include "integrate.h"
#include "rng.h"

typedef struct {
    u64 nK, nT;
    f64 S0;
    f64 *forward;         /* size = nT */
    f64 *maturities;      /* size = nT */
    u64 *maturities_days; /* size = nT */
    f64 *strikes;         /* size = nK */
    f64 *cprices;         /* size = nK x nT */
} pricing_surface;

typedef struct {
    f64 kappa; /* mean-reversion speed */
    f64 theta; /* long-run variance */
    f64 sigma; /* vol-of-vol */
    f64 rho;   /* correlation */
    f64 v0;    /* initial variance */
} heston_params;

////////////////////////////////////////////////////////////
//~ MC Pricing
////////////////////////////////////////////////////////////

internal void price_call_exact_sim(const heston_params *p, pricing_surface *s,
                                   Ran *rng, u64 M);

////////////////////////////////////////////////////////////
//~ Fourrier Pricing
////////////////////////////////////////////////////////////

internal f64c heston_cf(const heston_params *p, f64c u, f64 T);

////////////////////////////////////////////////////////////
//~ Lewis formula

internal f64 lewis_call_gl_quad(const heston_params *p, f64 F, f64 K, f64 T,
                                const gl_quad *q);

#endif
