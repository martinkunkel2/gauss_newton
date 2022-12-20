#include <gtest/gtest.h>
#include <gauss_newton.h>
#include <cmath>

struct ProblemWithJac : GaussNewton::Problem
{
    ProblemWithJac() : GaussNewton::Problem(2, 4)
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

struct ProblemWithoutJac : GaussNewton::Problem
{
    ProblemWithoutJac() : GaussNewton::Problem(2, 4)
    {
    }
    void eval(double *f, const double *x) const
    {
        f[0] = sin(x[0]);
        f[1] = sin(x[0] + 0.1);
        f[2] = sin(x[1]);
        f[3] = sin(x[1] - 0.1);
    }
};

TEST(GaussNewtonTest, TestJac)
{
    ProblemWithJac p1;
    ProblemWithoutJac p2;
    double x[2] = {0.2, 0.2};
    double J1[2 * 4];
    double J2[2 * 4];
    p1.jacobian(J1, x);
    p2.jacobian(J2, x);

    for (size_t i = 0; i < 8; i++)
        EXPECT_NEAR(J1[i], J2[i], 0.01);
}
