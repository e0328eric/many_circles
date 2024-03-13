#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_bessel.h>

#include "bessel_func.h"

#define I (gsl_complex_rect(0, 1))

/* FORTRAN functions in library */
/* Bessel function of the first kind. */
// NOTE: copied from zbsubs.for (the original source code)
//
// INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
//   ZR,ZI  - Z=CMPLX(ZR,ZI),  -PI.LT.ARG(Z).LE.PI
//   FNU    - ORDER OF INITIAL J FUNCTION, FNU.GE.0.0D0
//   KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
//            KODE= 1  RETURNS
//                     CY(I)=J(FNU+I-1,Z), I=1,...,N
//                = 2  RETURNS
//                     CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y)), I=1,...,N
//   N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
//
// OUTPUT     CYR,CYI ARE DOUBLE PRECISION
//   CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
//            CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
//            CY(I)=J(FNU+I-1,Z)  OR
//            CY(I)=J(FNU+I-1,Z)EXP(-ABS(Y))  I=1,...,N
//            DEPENDING ON KODE, Y=AIMAG(Z).
//   NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
//            NZ= 0   , NORMAL RETURN
//            NZ.GT.0 , LAST NZ COMPONENTS OF CY SET  ZERO DUE
//                      TO UNDERFLOW, CY(I)=CMPLX(0.0D0,0.0D0),
//                      I = N-NZ+1,...,N
//   IERR   - ERROR FLAG
//            IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
//            IERR=1, INPUT ERROR   - NO COMPUTATION
//            IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
//                    TOO LARGE ON KODE=1
//            IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
//                    BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
//                    REDUCTION PRODUCE LESS THAN HALF OF MACHINE
//                    ACCURACY
//            IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
//                    TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
//                    CANCE BY ARGUMENT REDUCTION
//            IERR=5, ERROR              - NO COMPUTATION,
//                    ALGORITHM TERMINATION CONDITION NOT MET
extern zbesj_wrap_t zbesj_wrap;

/* Bessel function of the second kind. */
// NOTE: copied from zbsubs.for (the original source code)
//
// INPUT      ZR,ZI,FNU ARE DOUBLE PRECISION
//   ZR,ZI  - Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
//            -PI.LT.ARG(Z).LE.PI
//   FNU    - ORDER OF INITIAL Y FUNCTION, FNU.GE.0.0D0
//   KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
//            KODE= 1  RETURNS
//                     CY(I)=Y(FNU+I-1,Z), I=1,...,N
//                = 2  RETURNS
//                     CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
//                     WHERE Y=AIMAG(Z)
//   N      - NUMBER OF MEMBERS OF THE SEQUENCE, N.GE.1
//   CWRKR, - DOUBLE PRECISION WORK VECTORS OF DIMENSION AT
//   CWRKI    AT LEAST N
//
// OUTPUT     CYR,CYI ARE DOUBLE PRECISION
//   CYR,CYI- DOUBLE PRECISION VECTORS WHOSE FIRST N COMPONENTS
//            CONTAIN REAL AND IMAGINARY PARTS FOR THE SEQUENCE
//            CY(I)=Y(FNU+I-1,Z)  OR
//            CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
//            DEPENDING ON KODE.
//   NZ     - NZ=0 , A NORMAL RETURN
//            NZ.GT.0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
//            UNDERFLOW (GENERALLY ON KODE=2)
//   IERR   - ERROR FLAG
//            IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
//            IERR=1, INPUT ERROR   - NO COMPUTATION
//            IERR=2, OVERFLOW      - NO COMPUTATION, FNU IS
//                    TOO LARGE OR CABS(Z) IS TOO SMALL OR BOTH
//            IERR=3, CABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
//                    BUT LOSSES OF SIGNIFCANCE BY ARGUMENT
//                    REDUCTION PRODUCE LESS THAN HALF OF MACHINE
//                    ACCURACY
//            IERR=4, CABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
//                    TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
//                    CANCE BY ARGUMENT REDUCTION
//            IERR=5, ERROR              - NO COMPUTATION,
//                    ALGORITHM TERMINATION CONDITION NOT MET
extern zbesy_wrap_t zbesy_wrap;

gsl_complex besselj(int nu, gsl_complex z) {
    // If nu < 0, then use the reflection formula
    if (nu < 0) {
        return gsl_complex_mul_real(besselj(-nu, z),
                                    (double)(1 - ((nu & 1) << 1)));
    }

    const double zr = GSL_REAL(z), zi = GSL_IMAG(z);
    const int kode = 1, n = 1;

    // output parameters
    double cyr, cyi;
    int nz, ierr;

    (*zbesj_wrap)(zr, zi, nu, kode, n, &cyr, &cyi, &nz, &ierr);

    if (ierr) {
        fprintf(stderr, "ERROR: cannot calculated besselj: errcode: %d\n",
                ierr);
        exit(EXIT_FAILURE);
    }

    return gsl_complex_rect(cyr, cyi);
}

gsl_complex bessely(int nu, gsl_complex z) {
    // If nu < 0, then use the reflection formula
    if (nu < 0) {
        return gsl_complex_mul_real(bessely(-nu, z),
                                    (double)(1 - ((nu & 1) << 1)));
    }

    const double zr = GSL_REAL(z), zi = GSL_IMAG(z);
    const int kode = 1, n = 1;

    // output parameters
    double cyr, cyi, cwrkr, cwrki;
    int nz, ierr;

    (*zbesy_wrap)(zr, zi, nu, kode, n, &cyr, &cyi, &nz, &cwrkr, &cwrki, &ierr);

    if (ierr) {
        fprintf(stderr, "ERROR: cannot calculated bessely: errcode: %d\n",
                ierr);
        exit(EXIT_FAILURE);
    }

    return gsl_complex_rect(cyr, cyi);
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
