#ifndef __linux__
#error "since this program uses <dlfcn.h>, only linux can run this program"
#endif

#include <assert.h>
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GSL_COMPLEX_LEGACY 1
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

#include "bessel_func.h"
#include "circle.h"
#include "single_circle_solv.h"
#include "utility.h"

#define OMEGA 5.0
#define N 11

#define DRAW_RADIUS 5.0
#define DRAW_DELTA 0.05

// N should be odd number
#if (N < 0 || !(N & 1))
#error "N should be odd positive number"
#endif

zbesj_wrap_t zbesj_wrap = NULL;
zbesy_wrap_t zbesy_wrap = NULL;

static void main_calculation(void);
static const gsl_vector_complex* make_eigenfunction_wave(
    const Circle* circle,
    double omega,
    const gsl_complex* datas,
    size_t datas_len
);

int main(void) {
    char* error;
    void* complex_bessel_lib =
        dlopen("./complex_bessel/build/lib/libcomplex_bessel.so", RTLD_LAZY);

    if (!complex_bessel_lib) {
        fprintf(stderr, "%s\n", dlerror());
        exit(EXIT_FAILURE);
    }

    // Clear any existing error
    dlerror();

    zbesj_wrap = (zbesj_wrap_t)dlsym(complex_bessel_lib, "zbesj_wrap");
    zbesy_wrap = (zbesy_wrap_t)dlsym(complex_bessel_lib, "zbesy_wrap");

    if ((error = dlerror()) != NULL) {
        fprintf(stderr, "%s\n", error);
        exit(EXIT_FAILURE);
    }

    // run main function
    main_calculation();

    dlclose(complex_bessel_lib);
}

static void main_calculation(void) {
    printf("Solve the Helmholtz equation where ω = %.5f\n", OMEGA);

    // Initialization for solving single domain problem
    Circle circle = {
        .center = CREAL(0),
        .radius = 1.0,
        .rho = CREAL(1.0),
        .kappa = CREAL(1.0),
    };
    assert(gsl_complex_abs(circle.rho) >= TOLERANCE &&
           "circle.rho should not be zero");

    printf("Solving that PDE on the single circle where N=%zu and\n",
           (size_t)N);
    dump_circle(stdout, &circle);

    gsl_complex data[N] = {
        CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(0.0),
        CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(5.0),
    };

    InitialCond_Single icond =
        initial_cond_single_init(&circle, OMEGA, data, N, make_eigenfunction_wave);
    solve_single_circle(&icond);

    FILE* data_real_file = fopen("datas_real.dat", "w+");
    FILE* data_imag_file = fopen("datas_imag.dat", "w+");
    FILE* data_abs_file = fopen("datas_abs.dat", "w+");

    fprintf(data_real_file, "#  x          y       real_val\n");
    fprintf(data_imag_file, "#  x          y       imag_val\n");
    fprintf(data_abs_file, "#  x          y        abs_val\n");

    gsl_complex coord, val;
    for (double r = DRAW_RADIUS; r > 0.0; r -= DRAW_DELTA) {
        for (double theta = 0.0; theta < 2 * M_PI; theta += 0.02) {
            coord = gsl_complex_polar(r, theta);
            val = get_solution_value(&icond, r, theta);

            fprintf(data_real_file, "%.5f    %.5f    %.5f\n", GSL_REAL(coord),
                    GSL_IMAG(coord), GSL_REAL(val));
            fprintf(data_imag_file, "%.5f    %.5f    %.5f\n", GSL_REAL(coord),
                    GSL_IMAG(coord), GSL_IMAG(val));
            fprintf(data_abs_file, "%.5f    %.5f    %.5f\n", GSL_REAL(coord),
                    GSL_IMAG(coord), gsl_complex_abs(val));
        }
    }

    fclose(data_abs_file);
    fclose(data_imag_file);
    fclose(data_real_file);
    initial_cond_single_deinit(&icond);
}

#define LOOP_FOR(i, len) \
    for (int i = ((int)len >> 1) - (int)len + 1; i <= (int)len >> 1; ++i)

#define IDX(i, len) (i + (len >> 1))

static const gsl_vector_complex* make_eigenfunction_wave(
    const Circle* circle,
    double omega,
    const gsl_complex* datas,
    size_t datas_len
) {
    assert(gsl_complex_abs(circle->center) < TOLERANCE && "`center` should be 0");
    assert(datas_len & 1 && "`datas_len` should be odd");
    assert(omega >= 0 && "`omega` should be positive");

    gsl_vector_complex* output = gsl_vector_complex_alloc(datas_len * 2);

    LOOP_FOR(i, datas_len) {
        size_t idx = IDX(i, datas_len);

        gsl_vector_complex_set(
            output, idx,
            gsl_complex_mul(datas[idx], besselj(i, CREAL(omega * circle->radius))));
        gsl_vector_complex_set(
            output, idx + datas_len,
            gsl_complex_mul_real(
                gsl_complex_mul(datas[idx], besselj_deriv(i, CREAL(omega * circle->radius))),
                omega));
    }

    return output;
}
