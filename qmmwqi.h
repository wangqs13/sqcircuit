#ifndef QMMWQI_H
#define QMMWQI_H

#include <vector>
#include <queue>

#include "matrix.h"
#include "orthogonal.h"
#include "almost.h"

template<class T>
std::vector<T> operator + (std::vector<T> a, const T& b)
{
    a.push_back(b);
    return a;
}

template<class V>
// V = float, double
class QuantumMealyMachineWithQuantumInput
{
public:
    typedef Matrix<V, Complex<V> > matrix;

    class Configuration
    {
    public:
        std::vector<int> pi, a;
        matrix rho1, rho2;
        Configuration (const std::vector<int> &pi, const std::vector<int> &a, const matrix &rho1, const matrix &rho2) : pi(pi), a(a), rho1(rho1), rho2(rho2)
        {
        }
        ~Configuration()
        {
        }
        friend std::ostream &operator << (std::ostream &out, const Configuration &config)
        {
            out << " input = ";
            for (auto x : config.pi) out << x << " ";
            out << std::endl;
            out << "output = ";
            for (auto x : config.a) out << x << " ";
            out << std::endl;
            //out << config.rho1;
            //out << config.rho2;
            return out;
        }
    };

    int dim_in, dim_s, gamma;
    matrix unitary;
    std::vector<matrix> measure;

    QuantumMealyMachineWithQuantumInput(int dim_in, int dim_s, const matrix &unitary, const std::vector<matrix> &measure) : dim_in(dim_in), dim_s(dim_s), unitary(unitary), measure(measure)
    {
        assert(unitary.is_dimension(dim_in*dim_s));
        gamma = measure.size();
        for (int i = 0; i < gamma; ++ i)
            assert(measure[i].is_dimension(dim_in*dim_s));
    }

    QuantumMealyMachineWithQuantumInput(const QuantumMealyMachineWithQuantumInput &m) : dim_in(m.dim_in), dim_s(m.dim_s), gamma(m.gamma), unitary(m.unitary), measure(m.measure)
    {
    }

    ~QuantumMealyMachineWithQuantumInput()
    {
    }

    matrix apply(const matrix &rho, const matrix &sigma, int a) const
    {
        assert(rho.is_dimension(dim_s));
        assert(sigma.is_dimension(dim_in));
        return (measure[a]*unitary*tensor_product(sigma, rho)*unitary.hermitian()*measure[a].hermitian()).trace(dim_s);
    }

    friend std::pair<bool, Configuration> check(const QuantumMealyMachineWithQuantumInput &m1, const matrix &rho1, const QuantumMealyMachineWithQuantumInput &m2, const matrix &rho2, const std::vector<matrix> &basis, AlmostEqual<V> equal)
    {
        assert(m1.dim_in == m2.dim_in && m1.gamma == m2.gamma);
        int dim_in = m1.dim_in;
        int gamma = m1.gamma;
        assert(rho1.is_dimension(m1.dim_s));
        assert(rho2.is_dimension(m2.dim_s));
        int dim = basis.size();
        for (int i = 0; i < dim; ++ i)
            assert(basis[i].is_dimension(dim_in));
        std::queue<Configuration> Q;
        Q.push(Configuration({}, {}, rho1, rho2));
        Orthogonal<V, matrix> ortho(equal);
        int cnt = 0;
        while (!Q.empty())
        {
            Configuration config = Q.front();
            Q.pop();
            matrix rho = direct_sum(config.rho1, -config.rho2);
            if (ortho.is_independent(rho))
            {
                //std::cout << ++ cnt << std::endl;
                //std::cout << config;
                if (!equal(config.rho1.trace().x, config.rho2.trace().x))
                    return std::pair<bool, Configuration>(false, config);
                ortho.add(rho);
                for (int sigma = 0; sigma < dim; ++ sigma)
                    for (int a = 0; a < gamma; ++ a)
                        Q.push(Configuration(config.pi+sigma, config.a+a, m1.apply(config.rho1, basis[sigma], a), m2.apply(config.rho2, basis[sigma], a)));
            }
        }
        return std::pair<bool, Configuration>(true, Configuration({}, {}, rho1, rho2));
    }
};

#endif // QMMWQI_H
