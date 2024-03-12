#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GSL_COMPLEX_LEGACY 1
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

#include "circle.h"
#include "single_circle_solv.h"
#include "utility.h"

#define OMEGA 5.0
#define N 11

#define DRAW_RADIUS 10.0
#define DRAW_DELTA 0.01

// N should be odd number
#if (N < 0 || !(N & 1))
#error "N should be odd positive number"
#endif

int main(void) {
    printf("Solve the Helmholtz equation where ω = %.5f\n", OMEGA);

    // Initialization for solving single domain problem
    Circle circle = {
        .center = CREAL(0),
        .radius = 1.0,
        .rho = -7.0,
        .kappa = 1.0,
    };
    assert(fabs(circle.rho) >= TOLERANCE && "circle.rho should not be zero");

    printf("Solving that PDE on the single circle where N=%zu and\n",
           (size_t)N);
    dump_circle(stdout, &circle);

    gsl_complex data[N] = {
        CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(0.0),  CREAL(0.0),
        CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(0.0), CREAL(5.0),
    };

    InitialCond_Single icond =
        initial_cond_single_init(&circle, OMEGA, data, N);
    solve_single_circle(&icond);

    FILE* data_real_file = fopen("datas_real.dat", "w+");
    FILE* data_imag_file = fopen("datas_imag.dat", "w+");
    FILE* data_abs_file = fopen("datas_abs.dat", "w+");
    FILE* circle_file = fopen("circle.dat", "w+");

    fprintf(circle_file, "#  x          y\n");
    gsl_complex circle_path;
    for (double theta = 0.0; theta < 2 * M_PI; theta += 0.01) {
        circle_path = gsl_complex_polar(circle.radius, theta);
        fprintf(circle_file, "%.5f    %.5f\n", GSL_REAL(circle_path),
                GSL_IMAG(circle_path));
    }

    fprintf(data_real_file, "#  x          y       real_val\n");
    fprintf(data_imag_file, "#  x          y       imag_val\n");

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

    fclose(circle_file);
    fclose(data_abs_file);
    fclose(data_imag_file);
    fclose(data_real_file);
    initial_cond_single_deinit(&icond);
}
