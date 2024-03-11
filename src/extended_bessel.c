#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "extended_bessel.h"

#define I (gsl_complex_rect(0, 1))

static gsl_complex besselj_purely_imag(int n, double x);
static gsl_complex bessely_purely_imag(int n, double x);

gsl_complex besselj(int n, gsl_complex z) {
    const double real = GSL_REAL(z), imag = GSL_IMAG(z);

    switch ((fabs(real) >= TOLERANCE) << 1 | (fabs(imag) >= TOLERANCE)) {
    case 0:
        return gsl_complex_rect(0, 0);
    case 1:
        return besselj_purely_imag(n, imag);
    case 2:
        return gsl_complex_rect(gsl_sf_bessel_Jn(n, real), 0);
    default:
        assert(0 && "z should be either real or puerly imaginary");
        break;
    }
}

gsl_complex bessely(int n, gsl_complex z) {
    const double real = GSL_REAL(z), imag = GSL_IMAG(z);

    switch ((fabs(real) >= TOLERANCE) << 1 | (fabs(imag) >= TOLERANCE)) {
    case 0:
        return gsl_complex_rect(0, 0);
    case 1:
        return bessely_purely_imag(n, imag);
    case 2:
        return gsl_complex_rect(gsl_sf_bessel_Yn(n, real), 0);
    default:
        assert(0 && "z should be either real or puerly imaginary");
        break;
    }
}

gsl_complex hankel1(int n, gsl_complex z) {
    return gsl_complex_add(besselj(n, z), gsl_complex_mul(bessely(n, z), I));
}

// derivatives of the first Bessel function and the Hankel function
// can be implemented with the recurrence relation 10.6.1
// reference: https://dlmf.nist.gov/10.6

gsl_complex besselj_deriv(int n, gsl_complex z) {
    return gsl_complex_div_real(
        gsl_complex_sub(besselj(n - 1, z), besselj(n + 1, z)), 2);
}

gsl_complex hankel1_deriv(int n, gsl_complex z) {
    return gsl_complex_div_real(
        gsl_complex_sub(hankel1(n - 1, z), hankel1(n + 1, z)), 2);
}

static gsl_complex besselj_purely_imag(int n, double x) {
    gsl_complex rotation;
    switch (n & 3) {
    case 0:
        GSL_SET_COMPLEX(&rotation, 1, 0);
        break;
    case 1:
        GSL_SET_COMPLEX(&rotation, 0, 1);
        break;
    case 2:
        GSL_SET_COMPLEX(&rotation, -1, 0);
        break;
    case 3:
        GSL_SET_COMPLEX(&rotation, 0, -1);
        break;
    default:
        assert(0 && "unreachable");
        break;
    }

    return gsl_complex_mul_real(rotation, gsl_sf_bessel_In(n, x));
}

static gsl_complex bessely_purely_imag(int n, double x) {
    gsl_complex rotation;
    switch (n & 3) {
    case 0:
        GSL_SET_COMPLEX(&rotation, 1, 0);
        break;
    case 1:
        GSL_SET_COMPLEX(&rotation, 0, 1);
        break;
    case 2:
        GSL_SET_COMPLEX(&rotation, -1, 0);
        break;
    case 3:
        GSL_SET_COMPLEX(&rotation, 0, -1);
        break;
    default:
        assert(0 && "unreachable");
        break;
    }

    gsl_complex output1, output2;
    GSL_SET_COMPLEX(&output1, gsl_sf_bessel_In(n, x), 0);
    GSL_SET_COMPLEX(&output2, M_2_PI * gsl_sf_bessel_Kn(n, x), 0);

    output1 = gsl_complex_mul(output1, rotation);
    output1 = gsl_complex_mul(output1, I);

    output2 = gsl_complex_div(output2, rotation);

    return gsl_complex_sub(output1, output2);
}
