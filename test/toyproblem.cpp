#include <gtest/gtest.h>
#include <gauss_newton.h>
#include <cmath>

struct ToyProblem : GaussNewton::Problem
{
    ToyProblem() : GaussNewton::Problem(1, 2)
    {
    }
    void eval(double *f, const double *x) const
    {
        f[0] = sin(x[0]);
        f[1] = sin(x[0] + 0.1);
    }
    void jacobian(double *J, const double *x) const
    {
        J[0] = cos(x[0]);
        J[1] = cos(x[0] + 0.1);
    }
};

TEST(GaussNewtonTest, ToyProblem)
{
    ToyProblem p;
    GaussNewton::Solver solver(p);
    double x0[1] = {0.2};
    double x[1];
    solver.solve(x, x0, 10, 0.001);
    printf("%lf\n", x[0]);
    EXPECT_NEAR(x[0], -0.05, 0.01);
}
