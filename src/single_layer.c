#include <assert.h>
#include <complex.h>
#include <math.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>

#include "extended_bessel.h"
#include "single_layer.h"

gsl_complex single_layer_potential(const Circle* domain, int n,
                                   gsl_complex wave_number, double r) {
    assert((r > 0) && "the value `r` should be positive");

    if (r > domain->radius) {
        return gsl_complex_mul_imag(
            gsl_complex_mul(
                besselj(n, gsl_complex_mul_real(wave_number, domain->radius)),
                hankel1(n, gsl_complex_mul_real(wave_number, r))),
            -M_PI_2 * domain->radius);
    } else {
        return gsl_complex_mul_imag(
            gsl_complex_mul(
                hankel1(n, gsl_complex_mul_real(wave_number, domain->radius)),
                besselj(n, gsl_complex_mul_real(wave_number, r))),
            -M_PI_2 * domain->radius);
    }
}

// k_star_helper is just calculates [J_nH1_n]'(z)
static inline gsl_complex k_star_helper(int n, gsl_complex val) {
    gsl_complex val1 = gsl_complex_mul(besselj_deriv(n, val), hankel1(n, val));
    gsl_complex val2 = gsl_complex_mul(besselj(n, val), hankel1_deriv(n, val));

    return gsl_complex_add(val1, val2);
}

gsl_complex k_star(const Circle* domain, int n, gsl_complex wave_number) {
    gsl_complex val =
        k_star_helper(n, gsl_complex_mul_real(wave_number, domain->radius));
    return gsl_complex_mul_imag(gsl_complex_mul(val, wave_number),
                                -M_PI_4 * domain->radius);
}
