#ifndef _HURRYPENG_ELECTRICFIELD_HPP
#define _HURRYPENG_ELECTRICFIELD_HPP

#include "Vector3D.hpp"
#include <list>
#include <functional>

namespace HurryPeng
{

class ElectricField
{
private:
    Vector3D uniformField;
    std::list<std::function<Vector3D(Vector3D)>> customOverlays;

public:
    ElectricField() : uniformField(), customOverlays() {}

    static Vector3D calcPointChargeField(const Vector3D & coordQ, const long double & q, const Vector3D & coordTgt)
    {
        const static long double K = 9E9;
        Vector3D deltaR = coordTgt - coordQ;
        Vector3D field;
        if (deltaR.norm() != 0.0) field = K * q / deltaR.norm() / deltaR.norm() * deltaR.unit();
        return field;
    }

    ElectricField & clear()
    {
        uniformField = Vector3D::ZERO_VECTOR;
        customOverlays.clear();
        return *this;
    }

    ElectricField & overlayUniformField(const Vector3D & newUniformField)
    {
        uniformField += newUniformField;
        return *this;
    }

    ElectricField & overlayPointChargeField(const Vector3D & coordQ, const long double & q)
    {
        customOverlays.push_back([coordQ, q](const Vector3D & coordTgt) -> Vector3D
            { return calcPointChargeField(coordQ, q, coordTgt); });
        return *this;
    }

    ElectricField & ovarlayCustomField(std::function<Vector3D(Vector3D)> fieldFunc)
    {
        customOverlays.push_back(fieldFunc);
        return *this;
    }

    Vector3D get(const Vector3D & coordTgt) const
    {
        Vector3D sum = uniformField;
        for (const auto & fieldFunc : customOverlays) sum += fieldFunc(coordTgt);
        return sum;
    }

    bool isEmpty() const
    {
        return uniformField == Vector3D(0, 0, 0) && customOverlays.empty();
    }
}; // class ElectricField

} // namespace HurryPeng

#endif
