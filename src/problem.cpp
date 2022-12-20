#include "problem.h"
#include <cmath>

void GaussNewton::Problem::jacobian(double *J, const double *x) const
{
    const double h = sqrt(__DBL_EPSILON__);
    // TODO: move to global place, to not allocate on each call
    double *f = new double[num_equations()];
    double *fh = new double[num_equations()];
    double *xh = new double[num_variables()];
    for (size_t i = 0; i < num_variables(); i++)
        xh[i] = x[i];

    eval(f, x);
    for (size_t i = 0; i < num_variables(); i++)
    {
        xh[i] += h;
        eval(fh, xh);

        for (size_t j = 0; j < num_equations(); j++)
            J[j * num_variables() + i] = (fh[j] - f[j]) / h;

        xh[i] = x[i];
    }
    delete xh;
    delete f;
    delete fh;
}
