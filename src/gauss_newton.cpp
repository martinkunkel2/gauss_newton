#include "gauss_newton.h"

#include <cstdlib>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

GaussNewton::Solver::Solver(const Problem &_p) : p(_p)
{
    dx = gsl_vector_alloc(p.num_variables());
    f = gsl_vector_alloc(p.num_equations());
    JTf = gsl_vector_alloc(p.num_variables());
    J = gsl_matrix_alloc(p.num_equations(), p.num_variables());
    JTJ = gsl_matrix_alloc(p.num_variables(), p.num_variables());
}

GaussNewton::Solver::~Solver()
{
    gsl_vector_free(dx);
    gsl_vector_free(f);
    gsl_vector_free(JTf);
    gsl_matrix_free(J);
    gsl_matrix_free(JTJ);
}

bool GaussNewton::Solver::solve(double *x, const double *x0, size_t maxIterations, double minStepSquare) const
{
    // Init
    for (size_t i = 0; i < p.num_variables(); i++)
    {
        x[i] = x0[i];
    }

    // Iteration
    for (size_t iter = 0; iter < maxIterations; iter++)
    {
        double stepsize = 0.0;
        bool result = step(x, stepsize);
        // some error occured
        if (!result)
            return false;
        // no further progress
        if (stepsize < minStepSquare)
            return true;
    }

    return false;
}

bool GaussNewton::Solver::step(double *x, double &stepsize) const
{
    // eval function and jacobian
    p.eval(f->data, x);
    p.jacobian(J->data, x);

    // calculate J^T * f
    if (gsl_blas_dgemv(CBLAS_TRANSPOSE_t::CblasTrans, 1.0, J, f, 0.0, JTf))
        return false;

    // TODO: Levenberg-Marquart regularisation
    // calculate J^T * J
    if (gsl_blas_dgemm(CBLAS_TRANSPOSE_t::CblasTrans, CBLAS_TRANSPOSE_t::CblasNoTrans, 1.0, J, J, 0.0, JTJ))
        return false;

    // cholesky decomp of J^T * J = L L^T
    if (gsl_linalg_cholesky_decomp1(JTJ))
        return false;

    // forward and backward solve  dx = (L L^T)^(-1) * (J^T f), dx -step
    if (gsl_linalg_cholesky_solve(JTJ, JTf, dx))
        return false;

    gsl_vector_view gx = gsl_vector_view_array(x, p.num_variables());
    // x = x - dx
    // TODO: Armijo search, Wolfe conditions
    if (gsl_vector_axpby(-1.0, dx, 1.0, &gx.vector))
        return false;

    if (gsl_blas_ddot(&gx.vector, &gx.vector, &stepsize))
        return false;

    return true;
}
