#ifndef COMPLEX_H
#define COMPLEX_H

#include <cmath>
#include <iostream>

template<class T> // T = int, float, double, ...
class Complex
{
public:
    T x, y;
    Complex(T x = 0, T y = 0) : x(x), y(y)
    {
    }
    ~Complex()
    {
    }

    friend std::istream &operator >> (std::istream &in, Complex &a)
    {
        return in >> a.x >> a.y;
    }

    friend std::ostream &operator << (std::ostream &out, const Complex &a)
    {
        out << a.x << " " << a.y;
        /*
        if (a.y == 0)
            out << a.x;
        else if (a.x == 0)
            out << a.y << "i";
        else if (a.y > 0)
            out << a.x << "+" << a.y << "i";
        else
            out << a.x << "-" << -a.y << "i";
            */
        return out;
    }

    T norm() const
    {
        return sqrt(norm2());
    }
    T norm2() const
    {
        return x*x+y*y;
    }

    Complex conjugate() const&
    {
        return Complex(x, -y);
    }
    Complex conjugate() &&
    {
        y = -y;
        return *this;
    }

    Complex &operator = (const Complex &a)
    {
        x = a.x;
        y = a.y;
        return *this;
    }

    Complex &operator += (const Complex &a)
    {
        x += a.x;
        y += a.y;
        return *this;
    }
    Complex &operator -= (const Complex &a)
    {
        x -= a.x;
        y -= a.y;
        return *this;
    }
    Complex &operator *= (const Complex &a)
    {
        T tx(x), ty(y);
        x = tx*a.x-ty*a.y;
        y = tx*a.y+ty*a.x;
        return *this;
    }
    Complex &operator /= (const Complex &a)
    {
        T tx(x), ty(y);
        T l(a.norm2());
        x = (tx*a.x+ty*a.y)/l;
        y = (-tx*a.y+ty*a.x)/l;
        return *this;
    }

    bool operator == (const Complex& a)
    {
        return x == a.x && y == a.y;
    }
    bool operator != (const Complex& a)
    {
        return !(*this == a);
    }

    friend Complex operator + (const Complex &a, const Complex &b)
    {
        return Complex(a.x+b.x, a.y+b.y);
    }
    friend Complex operator + (Complex &&a, const Complex &b)
    {
        return a += b;
    }
    friend Complex operator + (Complex &&a, Complex &&b)
    {
        return a += b;
    }
    friend Complex operator + (const Complex &a, Complex &&b)
    {
        return b += a;
    }

    friend Complex operator - (const Complex &a)
    {
        return Complex(-a.x, -a.y);
    }
    friend Complex operator - (Complex &&a)
    {
        a.x = -a.x;
        a.y = -a.y;
        return a;
    }
    friend Complex operator - (const Complex &a, const Complex &b)
    {
        return Complex(a.x-b.x, a.y-b.y);
    }
    friend Complex operator - (Complex &&a, const Complex &b)
    {
        return a -= b;
    }
    friend Complex operator - (Complex &&a, Complex &&b)
    {
        return a -= b;
    }
    friend Complex operator - (const Complex &a, Complex &&b)
    {
        b.x = a.x-b.x;
        b.y = a.y-b.y;
        return b;
    }

    friend Complex operator * (const Complex &a, const Complex &b)
    {
        return Complex(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x);
    }
    friend Complex operator * (Complex &&a, const Complex &b)
    {
        return a *= b;
    }
    friend Complex operator * (Complex &&a, Complex &&b)
    {
        return a *= b;
    }
    friend Complex operator * (const Complex &a, Complex &&b)
    {
        return b *= a;
    }

    friend Complex operator / (const Complex &a, const Complex &b)
    {
        T l(b.norm2());
        return Complex((a.x*b.x+a.y*b.y)/l, (-a.x*b.y+a.y*b.x)/l);
    }
    friend Complex operator / (Complex &&a, const Complex &b)
    {
        return a /= b;
    }
    friend Complex operator / (Complex &&a, Complex &&b)
    {
        return a /= b;
    }
    friend Complex operator / (const Complex &a, Complex &&b)
    {
        T tx(b.x), ty(b.y);
        T l(b.norm2());
        b.x = (a.x*tx+a.y*ty)/l;
        b.y = (-a.x*ty+a.y*tx)/l;
        return b;
    }

    friend Complex exp(const Complex &a)
    {
        return exp(a.x)*Complex(cos(a.y), sin(a.y));
    }
};

#endif // COMPLEX_H
