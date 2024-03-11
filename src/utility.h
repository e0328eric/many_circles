#ifndef UTILITY_H_
#define UTILITY_H_

#if defined(__GSL_COMPLEX_MATH_H__)
#define CREAL(val) (gsl_complex_rect(val, 0))
#define CIMAG(val) (gsl_complex_rect(0, val))
#endif

#if defined(GSL_REAL) && defined(GSL_IMAG)
#define CPX_FMT "%.2f + %.2fi"
#define CPX_ARG(val) GSL_REAL(val), GSL_IMAG(val)
#endif

#define TOLERANCE 1e-10

#endif // UTILITY_H_
