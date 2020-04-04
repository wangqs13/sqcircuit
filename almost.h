#ifndef ALMOST_H
#define ALMOST_H

template<class V>
V fabs(const V &a)
{
    return a >= 0 ? a : -a;
}

template<class V> // V = float, double, ...
class AlmostEqual
{
public:
    V eps;
    AlmostEqual(V eps) : eps(eps)
    {
    }
    AlmostEqual(const AlmostEqual &e) : eps(e.eps)
    {
    }
    ~AlmostEqual()
    {
    }

    bool operator () (const V &a, const V &b) const
    {
        return fabs(a-b) < eps;
        // this approximation may be not good enough.
    }
};

#endif // ALMOST_H
