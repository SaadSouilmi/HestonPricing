#ifndef HESTON_H
#define HESTON_H

#include "base.h"
#include "integrate.h"

typedef struct {
    f64 kappa; /* mean-reversion speed */
    f64 theta; /* long-run variance */
    f64 sigma; /* vol-of-vol */
    f64 rho;   /* correlation */
    f64 v0;    /* initial variance */
} heston_params;

internal f64c heston_cf(const heston_params *p, f64c u, f64 T);

////////////////////////////////////////////////////////////
//~ Lewis formula

internal f64 lewis_upper_bound(const heston_params *p, f64 T);
internal f64 lewis_call_gl_quad(const heston_params *p, f64 F, f64 K, f64 T,
                                const gl_quad *q);

#endif
