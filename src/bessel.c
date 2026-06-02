// SPDX-License-Identifier: MIT
#include <complex.h>
#include <math.h>

#include "bessel.h"

#define BESSEL_ZCUT 18.0

// sum_{k>=0} (z^2/4)^k / (k! Gamma(nu+k+1))
internal f64c
cbessel_i_series_sum(f64 nu, f64c z) {
    f64c z2 = z * z * 0.25;
    f64c term = 1.0 / tgamma(nu + 1.0);
    f64c sum = term;
    for (i32 k = 1; k < 200; k++) {
        term *= z2 / ((f64)k * (nu + (f64)k));
        sum += term;
        if (cabs(term) <= cabs(sum) * 1e-16)
            break;
    }
    return sum;
}

// sum_{k>=0} (-1)^k a_k / z^k, truncated at the smallest term (divergent)
internal f64c
cbessel_i_asymp_sum(f64 nu, f64c z) {
    f64 mu = 4.0 * nu * nu;
    f64c term = 1.0;
    f64c sum = 1.0;
    f64 prev = HUGE_VAL;
    for (i32 k = 1; k < 64; k++) {
        f64 m = 2.0 * (f64)k - 1.0;
        term *= -(mu - m * m) / (8.0 * (f64)k * z);
        f64 mag = cabs(term);
        if (mag > prev)
            break;
        prev = mag;
        sum += term;
        if (mag <= cabs(sum) * 1e-16)
            break;
    }
    return sum;
}

internal f64c
cbessel_i(f64 nu, f64c z) {
    if (cabs(z) < BESSEL_ZCUT)
        return cpow(z * 0.5, nu) * cbessel_i_series_sum(nu, z);
    return cexp(z) / csqrt(2.0 * PI * z) * cbessel_i_asymp_sum(nu, z);
}

internal f64c
cbessel_i_scaled(f64 nu, f64c z) {
    if (cabs(z) < BESSEL_ZCUT)
        return cexp(-z) * cpow(z * 0.5, nu) * cbessel_i_series_sum(nu, z);
    return cbessel_i_asymp_sum(nu, z) / csqrt(2.0 * PI * z);
}
