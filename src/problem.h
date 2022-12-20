#ifndef __GAUSS_NEWTON_PROBLEM_H__
#define __GAUSS_NEWTON_PROBLEM_H__

#include <cstdlib>

namespace GaussNewton
{
    struct Problem
    {
    public:
        Problem(size_t num_variables, size_t num_equations)
            : m_num_variables(num_variables),
              m_num_equations(num_equations)
        {
        }
        size_t num_variables() const { return m_num_variables; }
        size_t num_equations() const { return m_num_equations; }
        virtual void eval(double *f, const double *x) const = 0;
        virtual void jacobian(double *J, const double *x) const = 0;

    private:
        const size_t m_num_variables;
        const size_t m_num_equations;
    };

}

#endif
