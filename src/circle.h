#ifndef CIRCLE_H_
#define CIRCLE_H_

#include <stdio.h>

#include <gsl/gsl_complex.h>

typedef struct {
    gsl_complex center;
    double radius;
    double rho;
    double kappa;
} Circle;

void circle_set_center(Circle* domain, double x, double y);
gsl_complex wave_number(const Circle* domain, double omega);
void dump_circle(FILE* stream, const Circle* circle);

#endif // CIRCLE_H_
