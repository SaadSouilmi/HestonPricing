#ifndef INTEGRATE_H
#define INTEGRATE_H

#include "base.h"

typedef struct {
    f64 *x;
    f64 *w;
    u16 N;
} gl_quad;

internal void compute_gauss_legendre_roots_and_weights(gl_quad *q);

#endif
