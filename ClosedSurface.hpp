#ifndef _HURRYPENG_CLOSEDSURFACE_HPP
#define _HURRYPENG_CLOSEDSURFACE_HPP

#include "Vector3D.hpp"
#include <functional>
#include <list>
#include <vector>

namespace HurryPeng
{

struct ClosedSurface
{
    struct ParamSurface
    {
        std::function<Vector3D(const double &, const double &)> func;
        std::function<bool(const double &, const double &)> domain; // Domain function within a unit square

        std::list<Vector3D> surfacePoints(int precision = 200) const
        {
            std::list<Vector3D> points;
            double unit = 1 / (double)precision;
            for (int i = 0; i <= precision; i++) for (int j = 0; j <= precision; j++)
            {
                double u = i * unit, v = j * unit;
                if (domain(u, v)) points.push_back(func(u, v));
            }
            return points;
        }

        const static std::function<bool(const double &, const double &)> DOMAIN_FULL;

        static std::function<Vector3D(const double &, const double &)>
            generatePlate(const Vector3D & origin, const Vector3D & vecU, const Vector3D & vecV)
        {
            return [origin, vecU, vecV](const double & u, const double & v) -> Vector3D
                { return origin + u * vecU + v * vecV; };
        }
    };

    std::function<bool(const Vector3D &)> implicitSurface;
    std::vector<ParamSurface> paramSurfaces;

    bool covers(const Vector3D & coord) const { return implicitSurface(coord); }
    bool operator()(const Vector3D & coord) const { return covers(coord); }

    std::list<Vector3D> surfacePoints(int precision = 200) const
    {
        std::list<Vector3D> points;
        for (const ParamSurface paramSurface : paramSurfaces)
            points.splice(points.end(), paramSurface.surfacePoints(precision));
        return points;
    }

    static ClosedSurface generateSphere(const Vector3D & centre = {0, 0, 0}, const double & r = 0.01)
    {
        auto sphereImplicit = [r, centre](const Vector3D & coord) -> bool
            { return (coord - centre).norm() <= r; };

        auto sphereParam = [r, centre](const double & u, const double & v) -> Vector3D
        {
            double phi = u * 2 * PI, theta = v * PI;
            return centre + Vector3D{r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)};
        };

        ClosedSurface surface;
        surface.implicitSurface = sphereImplicit;
        surface.paramSurfaces.push_back(ClosedSurface::ParamSurface{sphereParam, ParamSurface::DOMAIN_FULL});
        return surface;
    }

    static ClosedSurface generateCube(const Vector3D & centre = {0, 0, 0}, const double & l = 0.01)
    {
        auto cubeImplicit = [centre, l](const Vector3D & coord) -> bool
        {
            return coord.x >= centre.x - l / 2 && coord.x < centre.x + l / 2
                && coord.y >= centre.y - l / 2 && coord.y < centre.y + l / 2
                && coord.z >= centre.z - l / 2 && coord.z < centre.z + l / 2;
        };

        Vector3D lowVertex = centre - Vector3D({l / 2, l / 2, l / 2});
        Vector3D highVertex = centre + Vector3D({l / 2, l / 2, l / 2});

        ClosedSurface surface;
        surface.implicitSurface = cubeImplicit;
        surface.paramSurfaces.push_back(ParamSurface{ParamSurface::generatePlate(lowVertex, {l, 0, 0}, {0, l, 0}), ParamSurface::DOMAIN_FULL});
        surface.paramSurfaces.push_back(ParamSurface{ParamSurface::generatePlate(lowVertex, {l, 0, 0}, {0, 0, l}), ParamSurface::DOMAIN_FULL});
        surface.paramSurfaces.push_back(ParamSurface{ParamSurface::generatePlate(lowVertex, {0, l, 0}, {0, 0, l}), ParamSurface::DOMAIN_FULL});
        surface.paramSurfaces.push_back(ParamSurface{ParamSurface::generatePlate(highVertex, {-l, 0, 0}, {0, -l, 0}), ParamSurface::DOMAIN_FULL});
        surface.paramSurfaces.push_back(ParamSurface{ParamSurface::generatePlate(highVertex, {-l, 0, 0}, {0, 0, -l}), ParamSurface::DOMAIN_FULL});
        surface.paramSurfaces.push_back(ParamSurface{ParamSurface::generatePlate(highVertex, {0, -l, 0}, {0, 0, -l}), ParamSurface::DOMAIN_FULL});

        return surface;
    }
};

const std::function<bool(const double &, const double &)>
    ClosedSurface::ParamSurface::DOMAIN_FULL = [](const double & u, const double & v)
    { return true; };

} // namespace HurryPeng

#endif
