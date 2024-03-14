#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_complex_double.h>

#include "bessel_func.h"
#include "many_circle_solv.h"
#include "single_layer.h"
#include "utility.h"

#define LOOP_FOR(i, len) \
    for (int i = ((int)len >> 1) - (int)len + 1; i <= (int)len >> 1; ++i)

#define IDX(i, len) (i + (len >> 1))

Circles_Array circles_array_alloc(size_t len) {
    Circles_Array output;

    output.circles = malloc(sizeof(Circle) * len);
    output.len = len;

    return output;
}

void circles_array_free(const Circles_Array* circle_arr) {
    free(circle_arr->circles);
}

const Circle* circles_array_get(const Circles_Array* circle_arr, size_t idx) {
    if (!circle_arr || idx >= circle_arr->len) {
        return NULL;
    }

    return &circle_arr->circles[idx];
}

bool circles_array_set(
    const Circles_Array* circle_arr,
    size_t idx,
    const Circle* val
) {
    if (!circle_arr || !val || idx >= circle_arr->len) {
        return false;
    }

    memcpy(&circle_arr->circles[idx], val, sizeof(Circle));
    return true;
}

////////////////////////////////////////////////////////////////////////////////
// static functions declarations
static void addition_theorem(
    gsl_vector_complex* output,
    const InitialCond_Multi* icond,
    size_t from,
    size_t into,
    int n,
    gsl_complex k
);

static void fill_top_left_matrix(
    gsl_matrix_complex* mat,
    const InitialCond_Multi* icond,
    size_t row,
    size_t column
);
static void fill_top_right_matrix(
    gsl_matrix_complex* mat,
    const InitialCond_Multi* icond,
    size_t row,
    size_t column
);
static void fill_bottom_left_matrix(
    gsl_matrix_complex* mat,
    const InitialCond_Multi* icond,
    size_t row,
    size_t column
);
static void fill_bottom_right_matrix(
    gsl_matrix_complex* mat,
    const InitialCond_Multi* icond,
    size_t row,
    size_t column
);


InitialCond_Multi initial_cond_multi_init(
    Circles_Array circles,
    double omega,
    gsl_complex* datas,
    size_t datas_len,
    Init_wave_fnt_t make_init_wave
) {
    InitialCond_Multi output;

    output.circles = circles;
    output.omega = omega;
    output.make_init_wave = make_init_wave;
    output.datas = malloc(sizeof(gsl_complex) * datas_len);
    output.datas_len = datas_len;
    memcpy(output.datas, datas, sizeof(gsl_complex) * datas_len);

    output.phi_outside = NULL;
    output.psi_inside = NULL;

    return output;
}

void initial_cond_multi_deinit(InitialCond_Multi* icond) {
    gsl_vector_complex_free(icond->phi_outside);
    gsl_vector_complex_free(icond->psi_inside);
    free(icond->datas);
}

gsl_complex get_solution_value_multi(
    const InitialCond_Multi* icond,
    double r,
    double theta
) {
    assert(0 && "UNIMPLEMENTED");
}

void solve_multi_circle(InitialCond_Multi* icond) {
    const size_t datas_len = icond->datas_len;
    const size_t mat_size = 2 * icond->datas_len * icond->circles.len;

    gsl_vector_complex* initial = gsl_vector_complex_alloc(mat_size);
    gsl_vector_complex* x= gsl_vector_complex_alloc(mat_size);
    gsl_matrix_complex* mat_a = gsl_matrix_complex_calloc(mat_size, mat_size);

    // TODO: fill initial conditions
    // for (size_t i = 0; i < icond->circles.len; ++i) {
    //     const Circle* circle = circles_array_get(&icond->circles, i);
    // }

    // filling fields of the matrix mat_a
    for (size_t row = 0; row < icond->circles.len; ++row) {
        for (size_t column = 0; column < icond->circles.len; ++column) {
            fill_top_left_matrix(mat_a, icond, row, column);
            fill_top_right_matrix(mat_a, icond, row, column);
            fill_bottom_left_matrix(mat_a, icond, row, column);
            fill_bottom_right_matrix(mat_a, icond, row, column);
        }
    }

    gsl_matrix_complex_free(mat_a);
    gsl_vector_complex_free(x);
    gsl_vector_complex_free(initial);
}

static void fill_top_left_matrix(
    gsl_matrix_complex* mat,
    const InitialCond_Multi* icond,
    size_t row,
    size_t column
) {
    if (row == column) {
        const Circle* circle = circles_array_get(&icond->circles, row);
        LOOP_FOR(k, icond->datas_len) {
            size_t idx = IDX(k, icond->datas_len) + row * icond->datas_len;

            gsl_matrix_complex_set(
                mat, idx, idx,
                gsl_complex_negative(single_layer_potential(
                    circle, k, CREAL(icond->omega), circle->radius)));
        }
    } else {
        // row: circle_into, column: circle_from because row is a row's idx and
        // column is a column's idx
        LOOP_FOR(k, icond->datas_len) {
            size_t idx = IDX(k, icond->datas_len);
            gsl_vector_complex_view vec_view = gsl_matrix_complex_subcolumn(
                    mat, icond->datas_len * column + k, icond->datas_len * row, icond->datas_len);
            addition_theorem(&vec_view.vector, icond, column, row, k, CREAL(icond->omega));
        }
    }
}

static void fill_top_right_matrix(
    gsl_matrix_complex* mat,
    const InitialCond_Multi* icond,
    size_t row,
    size_t column
) {
    const Circle* circle = circles_array_get(&icond->circles, row);
    gsl_complex wn = wave_number(circle, icond->omega);

    if (row == column) {
        LOOP_FOR(k, icond->datas_len) {
            size_t idx = IDX(k, icond->datas_len) + row * icond->datas_len;

            gsl_matrix_complex_set(
                mat, idx, idx + icond->datas_len * icond->circles.len,
                single_layer_potential(circle, k, wn, circle->radius));
        }
    } else {
        // row: circle_from, column: circle_into
        LOOP_FOR(k, icond->datas_len) {
            size_t idx = IDX(k, icond->datas_len);
            gsl_vector_complex_view vec_view = gsl_matrix_complex_subcolumn(
                mat, icond->datas_len * (column + icond->circles.len) + k,
                icond->datas_len * row, icond->datas_len);
            addition_theorem(&vec_view.vector, icond, row, column, k, wn);
        }
    }
}

static void fill_bottom_left_matrix(
    gsl_matrix_complex* mat,
    const InitialCond_Multi* icond,
    size_t row,
    size_t column
) {
    const Circle* circle = circles_array_get(&icond->circles, row);
    gsl_complex wn = wave_number(circle, icond->omega);

    if (row == column) {
        LOOP_FOR(k, icond->datas_len) {
            size_t idx = IDX(k, icond->datas_len) + row * icond->datas_len;

            gsl_matrix_complex_set(
                mat, idx + icond->datas_len * icond->circles.len, idx,
                gsl_complex_negative(gsl_complex_add_real(
                    k_star(circle, row, CREAL(icond->omega)), 0.5)));
        }
    } else {
        // row: circle_from, column: circle_into
        LOOP_FOR(k, icond->datas_len) {
            size_t idx = IDX(k, icond->datas_len);
            gsl_vector_complex_view vec_view = gsl_matrix_complex_subcolumn(
                mat, icond->datas_len * (column + icond->circles.len) + k,
                icond->datas_len * row, icond->datas_len);
            addition_theorem(&vec_view.vector, icond, row, column, k, wn);
        }
    }
}

static void addition_theorem(
    gsl_vector_complex* output,
    const InitialCond_Multi* icond,
    size_t from,
    size_t into,
    int n,
    gsl_complex wave_number
) {
    assert(output && "`output` cannot be NULL");
    assert(output->size == icond->datas_len && "length mismatched");

    const Circle *circle_from, *circle_into;
    if ((circle_from = circles_array_get(&icond->circles, from)) == NULL) {
        fprintf(stderr, "ERROR: index `from` out of bound\n");
        exit(1);
    }
    if ((circle_into = circles_array_get(&icond->circles, into)) == NULL) {
        fprintf(stderr, "ERROR: index `into` out of bound\n");
        exit(1);
    }

    gsl_complex centers_distance_vec =
        gsl_complex_sub(circle_into->center, circle_from->center);
    double centers_len = gsl_complex_abs(centers_distance_vec);
    double centers_arg = gsl_complex_arg(centers_distance_vec);

    // FIXME: I think this implementation is wrong!
    LOOP_FOR(l, icond->datas_len) {
        size_t idx = IDX(l, icond->datas_len);

        gsl_complex tmp1 = gsl_complex_mul(
            hankel1(n - l, gsl_complex_mul_real(wave_number, centers_len)),
            gsl_complex_polar(1, (double)n * centers_arg));

        gsl_complex tmp2 = gsl_complex_mul(
            besselj(l, gsl_complex_mul_real(wave_number, circle_from->radius)),
            gsl_complex_polar(1, (double)n * centers_arg));

        gsl_vector_complex_set(
            output, idx,
            gsl_complex_mul_real(gsl_complex_mul(tmp1, tmp2),
                1 - ((n - l & 1) >> 1))
        );
    }
}
