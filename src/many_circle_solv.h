#ifndef MANY_CIRCLE_SOLV_H_
#define MANY_CIRCLE_SOLV_H_

#include <stdbool.h>

#define GSL_COMPLEX_LEGACY 1
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector_complex_double.h>

#include "circle.h"
#include "initial_wave.h"

typedef struct {
    Circle* circles;
    size_t len;
} Circles_Array;

Circles_Array circles_array_alloc(size_t len);
void circles_array_free(const Circles_Array* circle_arr);

const Circle* circles_array_get(const Circles_Array* circle_arr, size_t idx);
bool circles_array_set(
    const Circles_Array* circle_arr,
    size_t idx,
    const Circle* val
);

typedef struct {
    Circles_Array circles;
    double omega;
    Init_wave_fnt_t make_init_wave;
    gsl_complex* datas;
    size_t datas_len;
    gsl_vector_complex* phi_outside;
    gsl_vector_complex* psi_inside;
} InitialCond_Multi;

InitialCond_Multi initial_cond_multi_init(
    Circles_Array circles,
    double omega,
    gsl_complex* datas,
    size_t datas_len,
    Init_wave_fnt_t make_init_wave
);
void initial_cond_multi_deinit(InitialCond_Multi* init_cond);

void solve_multi_circle(InitialCond_Multi* init_cond);
gsl_complex get_solution_value_multi(
    const InitialCond_Multi* init_cond,
    double r,
    double theta
);

#endif // MANY_CIRCLE_SOLV_H_
