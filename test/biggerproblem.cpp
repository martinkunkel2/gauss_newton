#include <gtest/gtest.h>
#include <gauss_newton.h>
#include <cmath>

struct BiggerProblem : GaussNewton::Problem
{
    BiggerProblem() : GaussNewton::Problem(11, 33)
    {
    }
    void eval(double *f, const double *x) const
    {
        for (int i = 0; i < 33; i++)
            f[i] = sin(x[i / 3] - i / 10.0);
    }
};

TEST(GaussNewtonTest, BiggerProblem)
{
    BiggerProblem p;
    GaussNewton::Solver solver(p);
    double x0[11] = {0.2};
    double x[11];
    solver.solve(x, x0, 10, 0.001);
}
