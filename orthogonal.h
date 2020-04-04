#ifndef ORTHOGONAL_H
#define ORTHOGONAL_H

#include <vector>

#include "almost.h"

template<class V, class T>
// e.g. V = double, T = Matrix<Complex<double> >
class Orthogonal
{
public:
    std::vector<T> basis;
    AlmostEqual<V> equal;
    Orthogonal(const AlmostEqual<V> &equal) : basis(), equal(equal)
    {
    }
    ~Orthogonal()
    {
    }

    class Gram_Schmidt
    {
    public:
        std::vector<T> basis;
        AlmostEqual<V> equal;
        Gram_Schmidt(const AlmostEqual<V>& equal) : basis(), equal(equal)
        {
        }
        ~Gram_Schmidt()
        {
        }
        void add(T a)
        {
            // Modified Gram-Schmidt Process.
            for (const auto& b : basis)
                a -= dot(b, a) * b;
            a.normalize();
            basis.push_back(std::move(a));
        }

        bool is_independent(const T& a)
        {
            V ret = 0;
            for (const auto& b : basis)
                ret += dot(b, a).norm2();
            return !equal(ret, a.norm2());
        }
    };

    void add(const T& a)
    {
        basis.push_back(a);
    }

    bool is_independent(const T& a)
    {
        Gram_Schmidt ortho(equal);
        for (const auto& b : basis)
            ortho.add(b);
        return ortho.is_independent(a);
        /*
        auto v = basis;
        int n = v.size();
        for (int i = 0; i < n; ++ i)
        {
            for (int j = 0; j < i-1; ++ j)
                v[i] -= dot(v[j], v[i])*v[j];
            v[i].normalize();
        }
        V ret = 0;
        for (const auto &b : v)
            ret += dot(b, a).norm2();
        std::cout << ret << " " << a.norm2() << std::endl;
        return !equal(ret, a.norm2());
        */
    }

};

#endif // ORTHOGONAL_H
