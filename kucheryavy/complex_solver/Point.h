#ifndef POINT_H_
#define POINT_H_

#include <iostream>
#include <cmath>
#include <string.h>

#define PI 3.14159265358979323846

namespace NR
{

class complex
{
    double real;
    double imag;
public:
    complex()
    {
        real = 0;
        imag = 0;
    }

    complex(double i_real, double i_imag)
    {
        real = i_real;
        imag = i_imag;
    }

    complex(const complex& n)
    {
        real = n.real;
        imag = n.imag;
    }

    double abs()
    {
        return sqrt(real * real + imag * imag);
    }

    complex operator+(const complex& n) const
    {
        complex sum;
        sum.real = this->real + n.real;
        sum.imag = this->imag + n.imag;
        return sum;
    }

    complex operator+(const double n) const
    {
        complex sum;
        sum.real = this->real + n;
        sum.imag = this->imag;
        return sum;
    }

    friend complex operator+(const double left, const complex& right)
    {
        complex sum;
        sum.real = left + right.real;
        sum.imag = right.imag;
        return sum;
    }

    complex operator-(const complex& n) const
    {
        complex diff;
        diff.real = this->real - n.real;
        diff.imag = this->imag - n.imag;
        return diff;
    }

    complex operator-(const double& n) const
    {
        complex diff;
        diff.real = this->real - n;
        diff.imag = this->imag;
        return diff;
    }

    complex operator*(const complex& n) const
    {
        complex mul;
        mul.real = this->real * n.real - this->imag * n.imag;
        mul.imag = this->real * n.imag + this->imag * n.real;
        return mul;
    }
    
    complex operator*(const double n) const
    {
        complex mul;
        mul.real = this->real * n;
        mul.imag = this->imag * n;
        return mul;
    }

    complex operator/(const double n) const
    {
        complex division;
        division.real = this->real / n;
        division.imag = this->imag / n;
        return division;
    }

    friend std::ostream& operator<<(std::ostream& stream, const complex n)
    {
        stream << n.real << " + " << n.imag << "i";
        return stream;
    }
};

/*
std::ostream &operator<<(std::ostream& stream, complex n)
{
    stream << n.real << " + " n.imag;
    return stream;
};*/

template <int dim>
class Point
{
    complex X[dim];
public:
    Point()
    {
        memset(X, 0x00, sizeof(complex)* dim);
    }

    Point(complex* i_X)
    {
        memcpy(X, i_X, sizeof(complex)* dim);
    }

    Point operator+(const Point& p) const
    {
        Point<dim>sum;
        for (int i = 0; i < dim; i++)
        {
            sum.X[i] = this->X[i] + p.X[i];
        }
        return sum;
    }

    Point operator+(double add) const
    {
        Point<dim>sum;
        for (int i = 0; i < dim; i++)
        {
            sum.X[i] = this->X[i] + add;
        }
        return sum;
    }

    Point& operator+=(const Point& p)
    {
        for (int i = 0; i < dim; i++)
        {
            this->X[i] += p.X[i];
        }
        return *this;
    }

    Point& operator+=(double add)
    {
        for (int i = 0; i < dim; i++)
        {
            this->X[i] += add;
        }
        return *this;
    }

    Point operator*(const double mul) const
    {
        Point<dim>res;
        for (int i = 0; i < dim; i++)
        {
            res.X[i] = this->X[i] * mul;
        }
        return res;
    }

    Point& operator*=(double mult)
    {
        for (int i = 0; i < dim; i++)
        {
            this->X[i] *= mult;
        }
        return *this;
    }

    Point operator-(Point &p)
    {
        Point<dim>sub;
        for (int i = 0; i < dim; i++)
        {
            sub.X[i] = this->X[i] - p.X[i];
        }
        return sub;
    }

    complex& operator[](int i)
    {
        return X[i];
    }

    const complex operator[](int i) const
    {
        return X[i];
    }

    double abs()
    {
        double sum = 0;
        for (int i = 0; i < dim; i++)
        {
            sum += X[i].abs() * X[i].abs();
        }
        return sqrt(sum);
    }


    template <int d>
    friend std::ostream &operator<<(std::ostream& stream, Point<d> p)
    {
        for (int i = 0; i < d; i++)
        {
            stream << p.X[i] << "\t";
        }
        return stream;
    };
};

}

#endif