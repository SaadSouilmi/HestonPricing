#include <complex.h>
#include <math.h>

#define LEWIS_INTEGRAL_TOL 1e-12

////////////////////////////////////////////////////////////
//~ MC Pricing
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
//~ Exact simulation

internal f64
real_integrated_var_cf(const heston_params *p, f64 V1, f64 V2, f64 dt, f64 d,
                       f64 z_kappa, f64c I_denom, f64 a) {
    f64c g = csqrt(p->kappa * p->kappa - 2 * p->sigma * p->sigma * I * a);
    f64 emkdt = exp(-p->kappa * dt);
    f64c emgdt = cexp(-g * dt);
    f64 svs = sqrt(V1 * V2);
    f64c z_gamma =
        svs * 4 * g * cexp(-0.5 * g * dt) / (p->sigma * p->sigma * (1 - emgdt));

    f64c term1 = g * cexp(-0.5 * (g - p->kappa) * dt) * (1 - emkdt) /
                 (p->kappa * (1 - emgdt));
    f64c term2 = cexp(
        (V1 + V2) / (p->sigma * p->sigma) *
        (p->kappa * (1 + emkdt) / (1 - emkdt) - g * (1 + emgdt) / (1 - emgdt)));
    f64c I_num = cbessel_i_scaled(0.5 * d - 1.0, z_gamma);
    f64c term3 = cexp(z_gamma - z_kappa) * I_num / I_denom;

    return creal(term1 * term2 * term3);
}

internal f64
sample_integrated_var(const heston_params *p) {}

internal void
price_call_exact_sim(const heston_params *p, pricing_surface *s, Ran *rng,
                     u64 M) {
    f64 V1 = p->v0;
    f64 V2 = 0;
    f64 S = p->S0;
    f64 t1 = 0;
    f64 t2 = 0;
    f64 dt = 0;
    for (u64 i = 0; i < s->nT; i++) {
        t2 = p->maturities[i];
        dt = t2 - t1;
        /* Sample V2 | V1 */
        f64 d = 4 * p->theta * p->kappa / (p->sigma * p->sigma);
        f64 emkdt = exp(-p->kappa * dt);
        f64 lambda =
            4 * p->kappa * emkdt / (p->sigma * p->sigma * (1 - emkdt)) * V1;
        V2 = ran_noncentral_chi2(rng, d, lambda);
        V2 *= (1 - emkdt) * p->sigma * p->sigma / (4 * p->kappa);
    }
}

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