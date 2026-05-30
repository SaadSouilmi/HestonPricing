#include <complex.h>
#include <math.h>

#define LEWIS_INTEGRAL_TOL 1e-12

////////////////////////////////////////////////////////////
//~ MC Pricing
////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////
//~ Fourrier Pricing
////////////////////////////////////////////////////////////

internal f64c
heston_cf(const heston_params *p, f64c u, f64 T) {
    f64c xi = p->kappa - I * p->rho * p->sigma * u;
    f64c d = csqrt(xi * xi + p->sigma * p->sigma * (I * u + u * u));
    f64c g = (xi - d) / (xi + d);

    f64c xi_minus_d = xi - d;
    f64c exp_mdT = cexp(-d * T);
    f64c one_minus_g_exp = 1 - g * exp_mdT;

    f64c C = p->kappa * p->theta / (p->sigma * p->sigma) *
             (xi_minus_d * T - 2 * clog(one_minus_g_exp / (1 - g)));
    f64c D =
        xi_minus_d * (1 - exp_mdT) / (p->sigma * p->sigma) / one_minus_g_exp;

    return cexp(C + D * p->v0);
}

////////////////////////////////////////////////////////////
//~ Lewis formula

internal f64
lewis_upper_bound(const heston_params *p, f64 T) {
    f64 A = sqrt(1.0 - p->rho * p->rho) / p->sigma *
            (p->kappa * p->theta * T + p->v0);
    f64 U = log(1.0 / LEWIS_INTEGRAL_TOL) / A;
    U = Max(U, 50.0);
    U = Min(U, 5000.0);
    return U;
}

internal f64
lewis_call_gl_quad(const heston_params *p, f64 F, f64 K, f64 T,
                   const gl_quad *q) {
    f64 U = lewis_upper_bound(p, T);
    f64 X = log(F / K);

    f64 integral = 0.0;
    for (u16 i = 0; i < q->N; i++) {
        f64c u = (q->x[i] + 1) * U / 2;
        f64c phi = heston_cf(p, u - 0.5 * I, T);
        integral +=
            creal(cexp(I * u * X) * phi / (u * u + 0.25)) * q->w[i] * U / 2;
    }

    return F - sqrt(F * K) / PI * integral;
}