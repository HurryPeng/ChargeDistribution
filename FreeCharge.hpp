#ifndef _HURRYPENG_FREECHARGE_HPP
#define _HURRYPENG_FREECHARGE_HPP

#include "Vector3D.hpp"
#include "ElectricField.hpp"

namespace HurryPeng
{

struct FreeCharge
{
    const static long double Q; // Unit charge that all free chargrs share
    const static long double M; // Unit mass that all free charges share

    bool isPositive;
    Vector3D coord;
    Vector3D vel = Vector3D::ZERO_VECTOR; // velocity

    FreeCharge(bool _isPositive, const Vector3D & _coord)
        :isPositive(_isPositive), coord(_coord) { }

    long double q() const { return isPositive ? Q : -Q; } // Returns signed q

    Vector3D accel(const ElectricField & field) const // acceleration
    {
        return q() * field.get(coord) / M;
    }

}; // struct FreeCharge

const long double FreeCharge::Q = 1.0E-12;
const long double FreeCharge::M = 1.0E-10;
    
} // namespace HurryPeng

#endif
