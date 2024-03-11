#ifndef SINGLE_CIRCLE_SOLV_H_
#define SINGLE_CIRCLE_SOLV_H_

#define GSL_COMPLEX_LEGACY 1
#include <gsl/gsl_complex.h>
#include <gsl/gsl_linalg.h>

#include "circle.h"

typedef struct {
    Circle circle;
    double omega;
    gsl_complex wave_number;
    gsl_complex* datas;
    size_t datas_len;
    gsl_vector_complex* phi_outside;
    gsl_vector_complex* psi_inside;
} InitialCond_Single;

InitialCond_Single initial_cond_single_init(const Circle* circle, double omega,
                                            gsl_complex* datas,
                                            size_t datas_len);
void initial_cond_single_deinit(InitialCond_Single* init_cond);

void solve_single_circle(InitialCond_Single* init_cond);
gsl_complex get_solution_value(const InitialCond_Single* init_cond, double r,
                               double theta);

#endif // SINGLE_CIRCLE_SOLV_H_
