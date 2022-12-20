#include <gtest/gtest.h>
#include <gauss_newton.h>
#include <cmath>

struct MultiVarProblem : GaussNewton::Problem
{
    MultiVarProblem() : GaussNewton::Problem(2, 4)
    {
    }
    void eval(double *f, const double *x) const
    {
        f[0] = sin(x[0]);
        f[1] = sin(x[0] + 0.1);
        f[2] = sin(x[1]);
        f[3] = sin(x[1] - 0.1);
    }
    void jacobian(double *J, const double *x) const
    {
        J[0] = cos(x[0]);
        J[1] = 0.0;
        J[2] = cos(x[0] + 0.1);
        J[3] = 0.0;
        J[4] = 0.0;
        J[5] = cos(x[1]);
        J[6] = 0.0;
        J[7] = cos(x[1] - 0.1);
    }
};

TEST(GaussNewtonTest, MultiVarProblem)
{
    MultiVarProblem p;
    GaussNewton::Solver solver(p);
    double x0[2] = {0.2, 0.2};
    double x[2];
    solver.solve(x, x0, 10, 0.001);
    printf("%lf\n", x[0]);
    printf("%lf\n", x[1]);
    EXPECT_NEAR(x[0], -0.05, 0.01);
    EXPECT_NEAR(x[1], 0.05, 0.01);
}
