#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

#include "extended_bessel.h"
#include "single_circle_solv.h"
#include "single_layer.h"
#include "utility.h"

#define LOOP_FOR(i, len) for (long i = (len >> 1) - len + 1; i <= len >> 1; ++i)

#define IDX(i, len) (i + (len >> 1))

static gsl_vector_complex* make_initial_wave(double r, double omega,
                                             gsl_complex* datas,
                                             size_t datas_len);

InitialCond_Single initial_cond_single_init(const Circle* circle, double omega,
                                            gsl_complex* datas,
                                            size_t datas_len) {
    assert(fabs(circle->rho) >= TOLERANCE && "circle.rho should not be zero");
    assert(omega >= 0 && "`omega` should be positive");

    InitialCond_Single output;

    memcpy(&output.circle, circle, sizeof(Circle));
    output.omega = omega;
    output.wave_number = wave_number(circle, omega);

    output.datas = malloc(sizeof(gsl_complex) * datas_len);
    output.datas_len = datas_len;
    memcpy(output.datas, datas, sizeof(gsl_complex) * datas_len);

    output.phi_outside = NULL;
    output.psi_inside = NULL;

    return output;
}

void initial_cond_single_deinit(InitialCond_Single* icond) {
    gsl_vector_complex_free(icond->psi_inside);
    gsl_vector_complex_free(icond->phi_outside);
    free(icond->datas);
}

void solve_single_circle(InitialCond_Single* icond) {
    const size_t n = icond->datas_len;
    const size_t mat_size = icond->datas_len * 2;

    gsl_vector_complex* initial = make_initial_wave(
        icond->circle.radius, icond->omega, icond->datas, icond->datas_len);
    gsl_matrix_complex* mat_a = gsl_matrix_complex_calloc(mat_size, mat_size);
    gsl_vector_complex* x = gsl_vector_complex_alloc(mat_size);

    LOOP_FOR(i, (long)n) {
        size_t idx = IDX(i, n);

        gsl_matrix_complex_set(
            mat_a, idx, idx,
            gsl_complex_negative(single_layer_potential(
                &icond->circle, i, CREAL(icond->omega), icond->circle.radius)));
        gsl_matrix_complex_set(mat_a, idx, idx + n,
                               single_layer_potential(&icond->circle, i,
                                                      icond->wave_number,
                                                      icond->circle.radius));
        gsl_matrix_complex_set(
            mat_a, idx + n, idx,
            gsl_complex_negative(gsl_complex_add_real(
                k_star(&icond->circle, i, CREAL(icond->omega)), 0.5)));
        gsl_matrix_complex_set(
            mat_a, idx + n, idx + n,
            gsl_complex_div_real(
                gsl_complex_add_real(
                    k_star(&icond->circle, i, icond->wave_number), -0.5),
                icond->circle.rho)); // assumed that ρ0 = 1
    }

    // using QR decomposition to solve Ax = init_wave.
    gsl_matrix_complex* mat_t = gsl_matrix_complex_alloc(mat_size, mat_size);
    gsl_linalg_complex_QR_decomp_r(mat_a, mat_t);
    gsl_linalg_complex_QR_solve_r(mat_a, mat_t, initial, x);

    gsl_vector_complex_const_view phi =
        gsl_vector_complex_const_subvector(x, 0, n);
    gsl_vector_complex_const_view psi =
        gsl_vector_complex_const_subvector(x, n, n);

    icond->phi_outside = gsl_vector_complex_alloc(n);
    icond->psi_inside = gsl_vector_complex_alloc(n);

    assert(icond->phi_outside && icond->psi_inside && "oom");

    gsl_vector_complex_memcpy(icond->phi_outside, &phi.vector);
    gsl_vector_complex_memcpy(icond->psi_inside, &psi.vector);

    gsl_vector_complex_free(x);
    gsl_matrix_complex_free(mat_t);
    gsl_matrix_complex_free(mat_a);
    gsl_vector_complex_free(initial);
}

gsl_complex get_solution_value(const InitialCond_Single* icond, double r,
                               double theta) {
    assert(icond->phi_outside && icond->psi_inside && "not solved yet");
    assert(r > 0 && "`r` should be positive");

    size_t n = icond->datas_len;
    gsl_complex output = GSL_COMPLEX_ZERO;

    if (r > icond->circle.radius) {
        LOOP_FOR(i, (long)n) {
            size_t idx = IDX(i, n);

            gsl_complex tmp =
                gsl_complex_mul(gsl_vector_complex_get(icond->phi_outside, idx),
                                single_layer_potential(&icond->circle, i,
                                                       CREAL(icond->omega), r));
            gsl_complex tmp2 = gsl_complex_mul(
                icond->datas[idx], besselj(i, CREAL(icond->omega * r)));

            tmp = gsl_complex_add(tmp, tmp2);
            tmp = gsl_complex_mul(tmp, gsl_complex_polar(1, i * theta));
            output = gsl_complex_add(output, tmp);
        }
    } else {
        LOOP_FOR(i, (long)n) {
            size_t idx = IDX(i, n);

            gsl_complex tmp =
                gsl_complex_mul(gsl_vector_complex_get(icond->psi_inside, idx),
                                single_layer_potential(&icond->circle, i,
                                                       icond->wave_number, r));
            tmp = gsl_complex_mul(tmp, gsl_complex_polar(1, i * theta));
            output = gsl_complex_add(output, tmp);
        }
    }

    return output;
}

static gsl_vector_complex* make_initial_wave(double r, double omega,
                                             gsl_complex* datas,
                                             size_t datas_len) {
    assert(datas_len & 1 && "`datas_len` should be odd");
    assert(omega >= 0 && "`omega` should be positive");

    gsl_vector_complex* output = gsl_vector_complex_alloc(datas_len * 2);

    LOOP_FOR(i, (long)datas_len) {
        size_t idx = IDX(i, datas_len);

        gsl_vector_complex_set(
            output, idx,
            gsl_complex_mul(datas[idx], besselj(i, CREAL(omega * r))));
        gsl_vector_complex_set(
            output, idx + datas_len,
            gsl_complex_mul_real(
                gsl_complex_mul(datas[idx], besselj_deriv(i, CREAL(omega * r))),
                omega));
    }

    return output;
}
