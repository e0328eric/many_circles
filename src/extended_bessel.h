#ifndef EXTENDED_BESSEL_H_
#define EXTENDED_BESSEL_H_

#define GSL_COMPLEX_LEGACY 1
#include <gsl/gsl_complex.h>

#ifndef TOLERANCE
#define TOLERANCE 1e-10
#endif

gsl_complex besselj(int n, gsl_complex z);
gsl_complex bessely(int n, gsl_complex z);
gsl_complex hankel1(int n, gsl_complex z);

gsl_complex besselj_deriv(int n, gsl_complex z);
gsl_complex hankel1_deriv(int n, gsl_complex z);

#endif // EXTENDED_BESSEL_H_
