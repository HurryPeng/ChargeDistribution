#ifndef _HURRYPENG_CONDUCTOR_HPP
#define _HURRYPENG_CONDUCTOR_HPP

#include "Vector3D.hpp"
#include "Chunk.hpp"
#include <functional>
#include <list>
#include <vector>
#include <map>
#include <tuple>
#include <float.h>
#include <random>

namespace HurryPeng
{

struct Conductor
{
    struct Surface
    {
        typedef std::function<Vector3D(const double &, const double &)> ParamSurface;
        typedef std::function<bool(const double &, const double &)> Domain;
        typedef std::function<bool(const Vector3D &)> ImplicitSurface;

        const static Domain DOMAIN_FULL;

        ParamSurface paramSurface;
        Domain domain; // Domain function within a unit square
        bool facingOut; // Whether the outer product for df/du and df/dv points out the surface

        ImplicitSurface implicitSurface;

        Vector3D at(const double & u, const double & v) const
        {
            return paramSurface(u, v);
        }
        Vector3D operator()(const double & u, const double & v) const
        {
            return this->at(u, v);
        }

        bool covers(const Vector3D & coord) const
        {
            return implicitSurface(coord);
        }
        bool operator()(const Vector3D & coord) const
        {
            return this->covers(coord);
        }

        std::list<std::pair<double, double>> validParams(int precision) const
        {
            double unit = 1.0 / precision;
            std::list<std::pair<double, double>> params;
            for (int i = 0; i < precision; i++)
                for (int j = 0; j < precision; j++)
                {
                    double u = i * unit, v = j * unit;
                    if (domain(u, v)) params.emplace_back(u, v);
                }
            return params;
        }

        std::list<Vector3D> surfacePoints(int precision) const
        {
            std::list<Vector3D> points;
            for (auto uv : validParams(precision))
                points.push_back(paramSurface(uv.first, uv.second));
            return points;
        }

        Vector3D normalVectorAt(const double & u, const double & v, int precision) const
        {
            double u1 = u + 1.0 / precision;
            double v1 = v + 1.0 / precision;

            Vector3D vecU = paramSurface(u1, v) - paramSurface(u, v);
            Vector3D vecV = paramSurface(u, v1) - paramSurface(u, v);

            Vector3D vecPerpendicular = vecU.outerProduct(vecV);
            if (!facingOut) vecPerpendicular = -vecPerpendicular;

            // Returns zero vector if normal vector cannot be calculated
            if (vecPerpendicular == Vector3D::ZERO_VECTOR) return Vector3D::ZERO_VECTOR;
            return vecPerpendicular.unit();
        }

        std::map<Vector3D, Vector3D> normalVectors(int precision) const
        {
            std::map<Vector3D, Vector3D> normVectors;
            for (auto uv : validParams(precision))
                normVectors[paramSurface(uv.first, uv.second)]
                    = normalVectorAt(uv.first, uv.second, precision);
            return normVectors;
        }

        static Surface generatePlate(const Vector3D & origin,
            const Vector3D & vecU, const Vector3D & vecV,
            bool facingOut = true, const Domain & domain = DOMAIN_FULL)
        {
            Surface plate;
            plate.paramSurface = [origin, vecU, vecV]
                (const double & u, const double & v) -> Vector3D
            {
                return origin + u * vecU + v * vecV;
            };
            plate.domain = domain;
            plate.facingOut = facingOut;
            plate.implicitSurface = [origin, vecU, vecV]
                (const Vector3D & coord) -> bool
            {
                return ((coord - origin) % (vecU * vecV) <= 0);
            };
            return plate;
        }
        
    };

    int precalcPrecision = 256;
    int spreadPrecision = 128;
    long double safetyDistance = 0.002;
    // Charges are spreaded this distance away from the actural surface
    // More usage later

    std::vector<Surface> surfaces;

    std::map<Chunk::ChunkId, std::map<int, Vector3D>> approxNormalVector;
    // approxNormalVector[chunkId][surfaceId]

    std::list<FreeCharge> boundCharges;

    std::vector<int> idOfSurfacesExceededBy(const Vector3D & coord) const
    {
        std::vector<int> ids;
        for (int i = 0; i < surfaces.size(); i++)
        {
            if (!surfaces[i].covers(coord)) ids.push_back(i);
        }
        return ids;
    }

    bool covers(const Vector3D & coord) const
    {
        return (idOfSurfacesExceededBy(coord).empty());
    }
    bool operator()(const Vector3D & coord) const { return covers(coord); }

    std::list<Vector3D> surfacePoints(int precision = 0) const
    {
        if (precision == 0) precision = precalcPrecision;
        std::list<Vector3D> points;
        for (const Surface & surface : surfaces)
            points.splice(points.end(), surface.surfacePoints(precision));
        return points;
    }

    void precalcApproxNormalVector(int regionRadius, int precision = 0)
    {
        if (precision == 0) precision = precalcPrecision;
        approxNormalVector.clear();

        for (int i = 0; i < surfaces.size(); i++)
        {
            const Surface & surface = surfaces[i];
            std::map<Vector3D, Vector3D> normalVectors = surface.normalVectors(precision);

            for (const Chunk::ChunkId & chunkId : Chunk::ChunkId::chunkIdsInRadius(regionRadius))
            {
                Vector3D centre = chunkId.centre();
                long double minDist = DBL_MAX;
                Vector3D normalVector;
                for (const std::pair<Vector3D, Vector3D> & pr : normalVectors)
                {
                    if (pr.second == Vector3D::ZERO_VECTOR) continue;
                    if ((pr.first - centre).norm() < minDist)
                    {
                        normalVector = pr.second;
                        minDist = (pr.first - centre).norm();
                    }
                }
                approxNormalVector[chunkId][i] = normalVector;
            }
        }
    }

    void spreadCharges(bool isPositive, int precision = 0)
    {
        if (precision == 0) precision = spreadPrecision;

        for (const Surface & surface : surfaces)
        {
            for (const std::pair<Vector3D, Vector3D> pr : surface.normalVectors(precision))
            {
                Vector3D surfacePoint = pr.first;
                Vector3D normalVector = -pr.second; // Pointing inwards
                if (normalVector == Vector3D::ZERO_VECTOR) continue;
                normalVector *= safetyDistance;
                surfacePoint += normalVector;
                if (!covers(surfacePoint)) continue;
                boundCharges.emplace_back(isPositive, surfacePoint);
            }
        }
    }

    std::list<std::pair<Vector3D, long double>> calcSurfaceField(int precision, long double distance) const
    {
        std::list<std::pair<Vector3D, long double>> fields;
        for (const Surface surface : surfaces)
        {
            for (const auto uv : surface.validParams(precision))
            {
                Vector3D surfacePoint = surface(uv.first, uv.second);
                Vector3D normalVect = surface.normalVectorAt(uv.first, uv.second, precalcPrecision);
                if (normalVect == Vector3D::ZERO_VECTOR) continue;
                Vector3D extendedPoint = surfacePoint + distance * normalVect;
                Vector3D intensity;
                for (const FreeCharge & freeCharge : boundCharges)
                {
                    Vector3D deltaX = extendedPoint - freeCharge.coord;
                    if (deltaX == Vector3D::ZERO_VECTOR) std::cerr << "?????\n";
                    intensity += freeCharge.q() / deltaX.norm() / deltaX.norm() * deltaX.unit();
                }
                fields.emplace_back(extendedPoint, intensity.norm());////////////////////
            }
        }
        return fields;
    }

    static Conductor generateEllipse
    (
        const Vector3D & centre = Vector3D::ZERO_VECTOR,
        const double & a = 0.01, const double & b = 0.01, const double & c = 0.01, 
        int precalcPrecision = 128, int spreadPrecision = 128
    )
    {
        Conductor conductor;
        Surface surface;
        
        surface.implicitSurface = [a, b, c, centre](const Vector3D & coord) -> bool
        {

            const auto & [x, y, z] = coord - centre;
            return x * x / a / a + y * y / b / b + z * z / c / c <= 1.0;
        };
        surface.paramSurface = [a, b, c, centre](const double & u, const double & v) -> Vector3D
        {
            double theta = v * PI, phi = u * 2 * PI;
            return centre + Vector3D{a * sin(theta) * cos(phi), b * sin(theta) * sin(phi), c * cos(theta)};
        };
        surface.domain = Surface::DOMAIN_FULL;
        surface.facingOut = false;

        conductor.surfaces.push_back(surface);
        conductor.precalcPrecision = precalcPrecision;
        conductor.spreadPrecision = spreadPrecision;

        return conductor;
    }

    static Conductor generateSphere
    (
        const Vector3D & centre = Vector3D::ZERO_VECTOR, const double & r = 0.01, 
        int precalcPrecision = 128, int spreadPrecision = 128
    )
    {
        /*
        Conductor conductor;
        Surface surface;
        surface.implicitSurface = [r, centre](const Vector3D & coord) -> bool
        {
            return (coord - centre).norm() <= r;
        };
        surface.paramSurface = [r, centre](const double & u, const double & v) -> Vector3D
        {
            double phi = u * PI, theta = v * 2 * PI;
            return centre + Vector3D{r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)};
        };
        surface.domain = Surface::DOMAIN_FULL;
        surface.facingOut = true;

        conductor.surfaces.push_back(surface);
        conductor.precalcPrecision = precalcPrecision;
        conductor.spreadPrecision = spreadPrecision;

        return conductor;
        */
        return generateEllipse(centre, r, r, r, precalcPrecision, spreadPrecision);
    }

    static Conductor generateCube(
        const Vector3D & centre = Vector3D::ZERO_VECTOR, const double & l = 0.01, 
        int precalcPrecision = 128, int spreadPrecision = 64)
    {
        Conductor conductor;
        Vector3D lowVertex = centre - Vector3D({l / 2, l / 2, l / 2});
        Vector3D highVertex = centre + Vector3D({l / 2, l / 2, l / 2});

        conductor.surfaces = 
        {
            Surface::generatePlate(lowVertex, {0, l, 0}, {l, 0, 0}), 
            Surface::generatePlate(lowVertex, {0, 0, l}, {0, l, 0}), 
            Surface::generatePlate(lowVertex, {l, 0, 0}, {0, 0, l}), 
            Surface::generatePlate(highVertex, {-l, 0, 0}, {0, -l, 0}), 
            Surface::generatePlate(highVertex, {0, -l, 0}, {0, 0, -l}),
            Surface::generatePlate(highVertex, {0, 0, -l}, {-l, 0, 0}) 
        };

        conductor.precalcPrecision = precalcPrecision;
        conductor.spreadPrecision = spreadPrecision;

        return conductor;
    }
};

const std::function<bool(const double &, const double &)>
    Conductor::Surface::DOMAIN_FULL = [](const double & u, const double & v)
    { return true; };

} // namespace HurryPeng

#endif
