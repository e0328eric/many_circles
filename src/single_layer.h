#ifndef SINGLE_LAYER_H_
#define SINGLE_LAYER_H_

#include <gsl/gsl_complex.h>

#include "circle.h"

gsl_complex single_layer_potential(const Circle* domain, int n,
                                   gsl_complex wave_number, double r);
gsl_complex k_star(const Circle* domain, int n, gsl_complex wave_number);

#endif // SINGLE_LAYER_H_
