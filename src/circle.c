#include <gsl/gsl_complex_math.h>

#include "circle.h"
#include "utility.h"

void circle_set_center(Circle* domain, double x, double y) {
    GSL_SET_COMPLEX(&domain->center, x, y);
}

gsl_complex wave_number(const Circle* domain, double omega) {
    return gsl_complex_sqrt_real(omega * omega * domain->rho / domain->kappa);
}

void dump_circle(FILE* stream, const Circle* circle) {
    fprintf(stream, "Circle {\n");
    fprintf(stream, "    " CPX_FMT "\n", CPX_ARG(circle->center));
    fprintf(stream, "    radius = %g\n", circle->radius);
    fprintf(stream, "    rho = %g\n", circle->rho);
    fprintf(stream, "    kappa = %g\n", circle->kappa);
    fprintf(stream, "}\n");
}