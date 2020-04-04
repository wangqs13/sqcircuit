#ifndef MATRIX_H
#define MATRIX_H

#include <initializer_list>
#include <algorithm>

#include "complex.h"

template<class V, class T>
// (V, T) = (float, float), (double, double)
//          (float, Complex<float>), (double, Complex<double>)
class Matrix
{
public:
    int row, col;
    T *value;
    T **a;
    Matrix(int row, int col, std::function<T (int, int)> p = [](int, int) -> T {return 0;}) : row(row), col(col)
    {
        assert(row > 0 && col > 0);
        value = new T[row*col]();
        a = new T*[row]();
        for (int i = 0; i < row; ++ i)
        {
            a[i] = value+i*col;
            for (int j = 0; j < col; ++ j)
                a[i][j] = p(i, j);
        }
    }
    Matrix(const std::initializer_list<std::initializer_list<T> > &v)
    {
        row = v.size();
        assert(row > 0);
        col = v.begin()->size();
        assert(col > 0);
        value = new T[row*col]();
        a = new T*[row]();
        int i = 0;
        for (auto &r : v)
        {
            assert(r.size() == col);
            a[i] = value+i*col;
            std::copy(r.begin(), r.end(), a[i]);
            ++ i;
        }
    }
    Matrix(const Matrix &A)
    {
        row = A.row;
        col = A.col;
        value = new T[row*col]();
        a = new T*[row]();
        for (int i = 0; i < row; ++ i)
            a[i] = value+i*col;
        std::copy(A.value, A.value+row*col, value);
    }
    Matrix(Matrix &&A)
    {
        row = A.row;
        col = A.col;
        value = A.value;
        a = A.a;
        A.value = NULL;
        A.a = NULL;
    }

    ~Matrix()
    {
        delete[] a;
        delete[] value;
    }

    friend std::ostream &operator << (std::ostream &out, const Matrix &a)
    {
        for (int i = 0; i < a.row; ++ i)
        {
            for (int j = 0; j < a.col; ++ j)
                out << a.a[i][j] << " ";
            out << std::endl;
        }
        return out;
    }

    T *operator [] (int index) const
    {
        return a[index];
    }

    Matrix transpose() const&
    {
        return Matrix(col, row, [&](int i, int j)
        {
            return a[j][i];
        });
    }
    Matrix transpose() &&
    {
        if (row == col)
        {
            for (int i = 0; i < row; ++ i)
                for (int j = i+1; j < row; ++ j)
                    std::swap(a[i][j], a[j][i]);
            return *this;
        }
        else
            return transpose();
    }

    Matrix conjugate() const&
    {
        return Matrix(row, col, [&](int i, int j)
        {
            return a[i][j].conjugate();
        });
        return *this;
    }
    Matrix conjugate() &&
    {
        for (int i = 0; i < row*col; ++ i)
            value[i] = value[i].conjugate();
        return *this;
    }

    Matrix hermitian() const&
    {
        return Matrix(col, row, [&](int i, int j)
        {
            return a[j][i].conjugate();
        });
    }
    Matrix hermitian() &&
    {
        if (row == col)
        {
            for (int i = 0; i < row; ++ i)
            {
                a[i][i] = a[i][i].conjugate();
                for (int j = i+1; j < row; ++ j)
                {
                    std::swap(a[i][j], a[j][i]);
                    a[i][j] = a[i][j].conjugate();
                    a[j][i] = a[j][i].conjugate();
                }
            }
            return *this;
        }
        else
            return hermitian();
    }

    bool is_dimension(int row, int col) const
    {
        return this->row == row && this->col == col;
    }
    bool is_dimension(int row) const
    {
        return this->row == row && this->col == row;
    }

    friend T dot(const Matrix &A, const Matrix &B)
    {
        assert(A.row == B.row && A.col == B.col);
        T ret = 0;
        for (int i = 0; i < A.row*A.col; ++ i)
            ret += A.value[i].conjugate()*B.value[i];
        return ret;
    }

    V norm2() const
    {
        V ret = 0;
        for (int i = 0; i < row*col; ++ i)
            ret += value[i].norm2();
        return ret;
    }
    V norm() const
    {
        return sqrt(norm2());
    }
    T trace() const
    {
        assert(row == col);
        T ret = 0;
        for (int i = 0; i < row; ++ i)
            ret += a[i][i];
        return ret;
    }
    Matrix trace(int dim) const
    {
        assert(row == col && row%dim == 0);
        return Matrix(dim, dim, [&](int x, int y)
        {
            T ret = 0;
            for (; x < row && y < row; x += dim, y += dim)
                ret += a[x][y];
            return ret;
        });
    }

    Matrix &normalize()
    {
        return *this /= norm();
    }

    Matrix &operator = (const Matrix &A)
    {
        if (this != &A)
        {
            delete[] a;
            delete[] value;
            row = A.row;
            col = A.col;
            value = new T[row * col]();
            a = new T * [row]();
            for (int i = 0; i < row; ++i)
                a[i] = value + i * col;
            std::copy(A.value, A.value + row * col, value);
        }
        return *this;
    }
    Matrix &operator = (Matrix &&A)
    {
        if (this != &A)
        {
            delete[] a;
            delete[] value;
            row = A.row;
            col = A.col;
            value = A.value;
            a = A.a;
            A.value = NULL;
            A.a = NULL;
        }
        return *this;
    }

    Matrix &operator += (const Matrix &A)
    {
        assert(row == A.row && col == A.col);
        for (int i = 0; i < row*col; ++ i)
            value[i] += A.value[i];
        return *this;
    }
    Matrix &operator -= (const Matrix &A)
    {
        assert(row == A.row && col == A.col);
        for (int i = 0; i < row*col; ++ i)
            value[i] -= A.value[i];
        return *this;
    }
    Matrix &operator *= (const Matrix &A)
    {
        assert(col == A.row);
        if (col == A.col)
        {
            T* tmp = new T[row*col];
            for (int i = 0; i < row; ++ i)
                for (int j = 0; j < col; ++ j)
                {
                    T &cur = tmp[i*col+j];
                    cur = 0;
                    for (int k = 0; k < col; ++ k)
                        cur += a[i][k]*A.a[k][j];
                }
            std::copy(tmp, tmp+row*col, value);
            delete[] tmp;
        }
        else
            *this = *this*A;
        return *this;
    }
    Matrix &operator *= (const T &number)
    {
        for (int i = 0; i < row*col; ++ i)
            value[i] *= number;
        return *this;
    }
    Matrix &operator /= (const T &number)
    {
        for (int i = 0; i < row*col; ++ i)
            value[i] /= number;
        return *this;
    }

    friend Matrix operator + (const Matrix &A, const Matrix &B)
    {
        assert(A.row == B.row && A.col == B.col);
        return Matrix(A.row, A.col, [&A, &B](int i, int j) -> T
        {
            return A.a[i][j]+B.a[i][j];
        });
    }
    friend Matrix operator + (Matrix &&A, const Matrix &B)
    {
        assert(A.row == B.row && A.col == B.col);
        return A += B;
    }
    friend Matrix operator + (Matrix &&A, Matrix &&B)
    {
        assert(A.row == B.row && A.col == B.col);
        return A += B;
    }
    friend Matrix operator + (const Matrix &A, Matrix &&B)
    {
        assert(A.row == B.row && A.col == B.col);
        return B += A;
    }

    friend Matrix operator - (const Matrix &A)
    {
        return Matrix(A.row, A.col, [&A](int i, int j) -> T
        {
            return -A.a[i][j];
        });
    }
    friend Matrix operator - (Matrix &&A)
    {
        for (int i = 0; i < A.row*A.col; ++ i)
            A.value[i] = -A.value[i];
        return A;
    }
    friend Matrix operator - (const Matrix &A, const Matrix &B)
    {
        assert(A.row == B.row && A.col == B.col);
        return Matrix(A.row, A.col, [&A, &B](int i, int j) -> T
        {
            return A.a[i][j]-B.a[i][j];
        });
    }
    friend Matrix operator - (Matrix &&A, const Matrix &B)
    {
        assert(A.row == B.row && A.col == B.col);
        return A -= B;
    }
    friend Matrix operator - (Matrix &&A, Matrix &&B)
    {
        assert(A.row == B.row && A.col == B.col);
        return A -= B;
    }
    friend Matrix operator - (const Matrix &A, Matrix &&B)
    {
        assert(A.row == B.row && A.col == B.col);
        return -(B -= A);
    }

    friend Matrix operator * (const Matrix &A, const Matrix &B)
    {
        assert(A.col == B.row);
        return Matrix(A.row, B.col, [&A, &B](int i, int j) -> T
        {
            T ret = 0;
            for (int k = 0; k < A.col; ++ k)
                ret += A.a[i][k]*B.a[k][j];
            return ret;
        });
    }
    friend Matrix operator * (Matrix &&A, const Matrix &B)
    {
        assert(A.col == B.row);
        if (A.col == B.col)
            return A *= B;
        else
            return A*B;
    }
    friend Matrix operator * (Matrix &&A, Matrix &&B)
    {
        assert(A.col == B.row);
        if (A.col == B.col)
            return A *= B;
        else
            return A*B;
    }
    friend Matrix operator * (const Matrix &A, Matrix &&B)
    {
        assert(A.col == B.row);
        if (A.row == B.row)
        {
            T* tmp = new T[A.row*B.col];
            for (int i = 0; i < A.row; ++ i)
                for (int j = 0; j < B.col; ++ j)
                {
                    T &cur = tmp[i*B.col+j];
                    cur = 0;
                    for (int k = 0; k < A.col; ++ k)
                        cur += A.a[i][k]*B.a[k][j];
                }
            std::copy(tmp, tmp+A.row*B.col, B.value);
            delete[] tmp;
            return B;
        }
        else
            return A*B;
    }

    friend Matrix operator * (const Matrix &A, const T &number)
    {
        return Matrix(A.row, A.col, [&A, &number](int i, int j) -> T
        {
            return number*A.a[i][j];
        });
    }
    friend Matrix operator * (Matrix &&A, const T &number)
    {
        return A *= number;
    }
    friend Matrix operator * (const T &number, const Matrix &A)
    {
        return A*number;
    }
    friend Matrix operator * (const T &number, Matrix &&A)
    {
        return A *= number;
    }

    friend Matrix operator / (const Matrix &A, const T &number)
    {
        return Matrix(A.row, A.col, [&A, &number](int i, int j) -> T
        {
            return A.a[i][j]/number;
        });
    }
    friend Matrix operator / (Matrix &&A, const T &number)
    {
        return A /= number;
    }

    friend Matrix direct_sum(const Matrix &A, const Matrix &B)
    {
        return Matrix(A.row+B.row, A.col+B.col, [&A, &B](int i, int j) -> T
        {
            if (i < A.row && j < A.col)
                return A.a[i][j];
            else if (i >= A.row && j >= A.col)
                return B.a[i-A.row][j-A.col];
            else
                return 0;
        });
    }

    bool operator == (const Matrix& A) const
    {
        if (!(row == A.row && col == A.col)) return false;
        for (int i = 0; i < row; ++i)
            for (int j = 0; j < col; ++j)
                if (a[i][j] != A.a[i][j]) return false;
        return true;
    }

    // this method is not efficient enough.
    template<class TA, class... TB>
    friend Matrix direct_sum(const TA &A, const TB&... B)
    {
        return direct_sum(A, direct_sum(B...));
    }

    friend Matrix tensor_product(const Matrix &A, const Matrix &B)
    {
        return Matrix(A.row*B.row, A.col*B.col, [&A, &B](int i, int j) -> T
        {
            int ia = i/B.row;
            int ib = i%B.row;
            int ja = j/B.col;
            int jb = j%B.col;
            return A.a[ia][ja]*B.a[ib][jb];
        });
    }
    template<class TA, class... TB>
    friend Matrix tensor_product(const TA &A, const TB&... B)
    {
        return tensor_product(A, tensor_product(B...));
    }

    static Matrix zeros(int row, int col)
    {
        return Matrix(row, col);
    }
    static Matrix zeros(int row)
    {
        return Matrix(row, row);
    }

    static Matrix eye(int row)
    {
        return Matrix(row, row, [](int i, int j) -> T {return i == j;});
    }

    static Matrix CNOT(int dim, int ctrl, int targ)
    {
        assert(ctrl != targ);
        assert(1 <= ctrl && ctrl <= dim);
        assert(1 <= targ && targ <= dim);
        return Matrix(1<<dim, 1<<dim, [&](int i, int j) -> T
        {
            j ^= (j>>(dim-ctrl)&1) << (dim-targ);
            return i == j;
        });
    }

    static Matrix Toffoli(int dim, int ctrl1, int ctrl2, int targ)
    {
        assert(ctrl1 != ctrl2 && targ != ctrl1 && targ != ctrl2);
        assert(1 <= ctrl1 && ctrl1 <= dim);
        assert(1 <= ctrl2 && ctrl2 <= dim);
        assert(1 <= targ && targ <= dim);
        return Matrix(1<<dim, 1<<dim, [&](int i, int j) -> T
        {
            j ^= ((j>>(dim-ctrl1)&1) && (j>>(dim-ctrl2)&1)) << (dim-targ);
            return i == j;
        });
    }

};

#endif // MATRIX_H
