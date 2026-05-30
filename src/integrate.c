#include <math.h>

#define GL_NEWTON_ITER 6
#define GL_NEWTON_TOL 1e-15

#define ComputeLegendrePoly(x, p0, p1, p2, N)                                  \
    do {                                                                       \
        p0 = 1;                                                                \
        p1 = x;                                                                \
        for (u16 j = 1; j < N; j++) {                                          \
            p2 = ((2 * j + 1) * x * p1 - j * p0) / (j + 1);                    \
            p0 = p1;                                                           \
            p1 = p2;                                                           \
        }                                                                      \
    } while (0)

internal void
compute_gauss_legendre_roots_and_weights(gl_quad *q) {
    u16 m = q->N / 2;
    f64 x0, x1, p0, p1, p2, p2_d;
    for (u16 i = 0; i < m; i++) {
        x0 = cos(PI * (i + 0.75) / (q->N + 0.5));
        for (u16 iter = 0; iter < GL_NEWTON_ITER; iter++) {
            ComputeLegendrePoly(x0, p0, p1, p2, q->N);
            p2_d = q->N * (p0 - x0 * p2) / (1 - x0 * x0);
            x1 = x0 - p2 / p2_d;
            if (fabs(x1 - x0) <= GL_NEWTON_TOL)
                break;
            x0 = x1;
        }
        ComputeLegendrePoly(x1, p0, p1, p2, q->N);
        p2_d = q->N * (p0 - x1 * p2) / (1 - x1 * x1);
        f64 w = 2 / (1 - x1 * x1) / (p2_d * p2_d);

        q->x[i] = -x1;
        q->x[q->N - 1 - i] = x1;
        q->w[i] = w;
        q->w[q->N - 1 - i] = w;
    }

    if (q->N & 1) {
        q->x[m] = 0;
        ComputeLegendrePoly(0, p0, p1, p2, q->N);
        p2_d = q->N * p0;
        q->w[m] = 2 / (p2_d * p2_d);
    }
}