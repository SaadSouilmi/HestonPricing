// SPDX-License-Identifier: MIT
#ifndef BESSEL_H
#define BESSEL_H

#include "base.h"

// I_nu(z), real order nu > -1, complex z (right half-plane).
internal f64c cbessel_i(f64 nu, f64c z);

// I_nu(z) * exp(-z), for overflow-safe ratios.
internal f64c cbessel_i_scaled(f64 nu, f64c z);

#endif
