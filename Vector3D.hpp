#ifndef _HURRYPENG_VECTOR3D_HPP
#define _HURRYPENG_VECTOR3D_HPP

#include <exception>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

namespace HurryPeng
{

const long double PI = 3.1415926535897932384626433832795;

struct Vector3D
{
    Vector3D() : x(0.0), y(0.0), z(0.0) {}
    Vector3D(const long double & _x, const long double & _y, const long double & _z)
        :x(_x), y(_y), z(_z) {}
    Vector3D(const std::initializer_list<long double> & il)
    {
        if (il.size() != 3) throw std::invalid_argument("Initializing a Vector3D with not exactly 3 arguments");
        x = il.begin()[0], y = il.begin()[1], z = il.begin()[2];
    }

    long double x;
    long double y;
    long double z;

    bool operator==(const Vector3D & rhs) const
    {
        return (x == rhs.x && y == rhs.y && z == rhs.z);
    }

    bool operator!=(const Vector3D & rhs) const
    {
        return !(*this == rhs);
    }

    Vector3D operator+(const Vector3D & rhs) const
    {
        return Vector3D(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    Vector3D operator-() const
    {
        return Vector3D(-x, -y, -z);
    }

    Vector3D operator-(const Vector3D & rhs) const
    {
        return *this + -rhs;
    }

    friend Vector3D operator*(const long double & lhs, const Vector3D & rhs);
    friend Vector3D operator*(const Vector3D & lhs, const long double & rhs);
    friend Vector3D operator/(const Vector3D & lhs, const long double & rhs);

    Vector3D & operator+=(const Vector3D & rhs)
    {
        return *this = *this + rhs;
    }

    Vector3D & operator-=(const Vector3D & rhs)
    {
        return *this = *this - rhs;
    }

    Vector3D & operator*=(const long double & rhs)
    {
        return *this = *this * rhs;
    }

    Vector3D & operator/=(const long double & rhs)
    {
        return *this = *this / rhs;
    }

    long double innerProduct(const Vector3D & rhs) const
    {
        return x * rhs.x + y * rhs.y + z * rhs.z;
    }

    long double operator%(const Vector3D & rhs) const
    {
        return innerProduct(rhs);
    }

    Vector3D outerProduct(const Vector3D & rhs) const
    {
        return Vector3D(y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x);
    }

    Vector3D operator*(const Vector3D & rhs) const
    {
        return outerProduct(rhs);
    }

    bool operator<(const Vector3D & rhs) const
    {
        if (x != rhs.x) return x < rhs.x;
        if (y != rhs.y) return y < rhs.y;
        if (z != rhs.z) return z < rhs.z;
        return false;
    }

    long double norm() const
    {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vector3D unit() const
    {
        if (this->norm() == 0.0) return ZERO_VECTOR; //throw std::invalid_argument("Calculating unit of a zero Vector3D");
        return *this / this->norm();
    }

    long double cosineOfAngle(const Vector3D & rhs) const
    {
        if (this->norm() == 0.0 || rhs.norm() == 0.0)
            throw std::invalid_argument("Calculating angle with a zero Vector3D");
        return *this % rhs / this->norm() / rhs.norm();
    }

    Vector3D projectOnto(const Vector3D & v) const
    {
        return this->innerProduct(v) * v.unit();
    }

    operator std::string() const
    {
        std::stringstream ss;
        ss << "{" << x << ", " << y << ", " << z << "}";
        return ss.str();
    }

    friend std::ostream & operator<<(std::ostream & lhs, const Vector3D & rhs);

    static const Vector3D ZERO_VECTOR;
}; // class Vector3D

Vector3D operator*(const long double & lhs, const Vector3D & rhs)
{
    return Vector3D(rhs.x * lhs, rhs.y * lhs, rhs.z * lhs);
}

Vector3D operator*(const Vector3D & lhs, const long double & rhs)
{
    return rhs * lhs;
}

Vector3D operator/(const Vector3D & lhs, const long double & rhs)
{
    if (rhs == 0.0) throw std::invalid_argument("Vector3D divided by zero");
    return Vector3D(lhs.x / rhs, lhs.y / rhs, lhs.z / rhs);
}

std::ostream & operator<<(std::ostream & lhs, const Vector3D & rhs)
{
    lhs << std::string(rhs);
    return lhs;
}

const Vector3D Vector3D::ZERO_VECTOR = {0, 0, 0};

} // namespace HurryPeng

#endif
