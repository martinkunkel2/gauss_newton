#ifndef __GAUSS_NEWTON_H__
#define __GAUSS_NEWTON_H__

#include "problem.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

namespace GaussNewton
{
    class Solver
    {
    public:
        Solver(const Problem &_p);
        ~Solver();
        bool solve(double *x, const double *x0, size_t maxIterations, double minStepSquare) const;

    private:
        bool step(double *x, double &stepobj, double &obj) const;
        const Problem &p;
        gsl_vector *dx;
        gsl_vector *f;
        gsl_vector *JTf;
        gsl_matrix *J;
        gsl_matrix *JTJ;
    };
}

#endif
