#ifndef INITIAL_WAVE_H_
#define INITIAL_WAVE_H_

#include <gsl/gsl_vector_complex_double.h>

#include "circle.h"

// NOTE:
// 1st param: const Circle* circle     (circle data)
// 2nd param: double omega             (wave frequency)
// 3rd param: const gsl_complex* datas (see the last note)
// 4th param: size_t datas_len         (datas pointer length)
//
// NOTE: this function should return length (datas_len * 2) amount of initial
// datas where first top (datas_len) part contains values of the initial wave
// value at the boundary of a circle, and the last bottom (datas_len) part
// contans values of its normal derivative.
//
// NOTE:
// datas contains constants c_j where u^in = 󰒠 c_j J_j(omega * r) e^(i j θ)
// where this representation is expressed with respect to the polar coordinate
// whose origin is 0. In other cases, please use the addition theorem to modify
// the datas before putting it into this function.
// In addition, that sum ranges over -(datas_len / 2) <= j <= (datas_len / 2).
typedef const gsl_vector_complex* (*Init_wave_fnt_t)(const Circle*, double, const gsl_complex*, size_t);

#endif // INITIAL_WAVE_H_
