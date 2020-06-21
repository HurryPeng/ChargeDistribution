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
#include <numeric>
#include <algorithm>
#include <optional>

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
        const static Domain DOMAIN_CIRCLE;

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
            for (auto [u, v] : validParams(precision))
                points.push_back(paramSurface(u, v));
            return points;
        }

        Vector3D rawNormalVectorAt(const double & u, const double & v, int precision) const
        {
            double u1 = u + 1.0 / precision;
            double v1 = v + 1.0 / precision;

            Vector3D du = paramSurface(u1, v) - paramSurface(u, v);
            Vector3D dv = paramSurface(u, v1) - paramSurface(u, v);

            Vector3D vecPerpendicular = du.outerProduct(dv);
            if (!facingOut) vecPerpendicular = -vecPerpendicular;

            // Returns zero vector if normal vector cannot be calculated
            if (vecPerpendicular == Vector3D::ZERO_VECTOR) return Vector3D::ZERO_VECTOR;
            return vecPerpendicular;
        }

        Vector3D normalVectorAt(const double & u, const double & v, int precision) const
        {
            return rawNormalVectorAt(u, v, precision).unit();
        }

        std::map<Vector3D, Vector3D> normalVectors(int precision) const
        {
            std::map<Vector3D, Vector3D> normVectors;
            for (auto [u, v] : validParams(precision))
                normVectors[paramSurface(u, v)]
                    = normalVectorAt(u, v, precision);
            return normVectors;
        }

        long double gaussianCurvatureAt(const double & u, const double & v, int precision) const
        {
            precision *= 16; // Very high precision is required because of the use of second-order derivatives
            double u1 = u + 0.5 / precision, u2 = u + 1.0 / precision;
            double v1 = v + 0.5 / precision, v2 = v + 1.0 / precision;

            Vector3D 
                du = paramSurface(u1, v) - paramSurface(u, v),
                du1 = paramSurface(u2, v) - paramSurface(u1, v),
                ddu = du1 - du,
                dv = paramSurface(u, v1) - paramSurface(u, v),
                dv1 = paramSurface(u, v2) - paramSurface(u, v1),
                ddv = dv1 - dv,
                du10 = paramSurface(u1, v1) - paramSurface(u, v1),
                du11 = paramSurface(u2, v1) - paramSurface(u1, v1),
                ddu1 = du11 - du10,
                dduv = ddu1 - ddu;
            Vector3D n = normalVectorAt(u, v, precision * 2);

            long double
                E = du % du, 
                F = du % dv,
                G = dv % dv,
                L = ddu % n,
                M = dduv % n,
                N = ddv % n;

            long double K = (L * N - M * M) / (E * G - F * F);
            return K;
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
            plate.implicitSurface = [origin, vecU, vecV, facingOut]
                (const Vector3D & coord) -> bool
            {
                if (facingOut) return ((coord - origin) % (vecU * vecV) <= 0);
                return ((coord - origin) % (vecU * vecV) >= 0);
            };
            return plate;
        }
    }; // struct Surface

    int precalcPrecision = 256;
    int spreadPrecision = 128;
    long double safetyDistance = 0.002;
    // Charges are spreaded this distance away from the actural surface
    // More usage later

    std::vector<Surface> surfaces;

    std::map<Chunk::ChunkId, std::map<int, std::pair<double, double>>> approxClosestSurfacePoint;
    // approxClosestSurfacePoint[chunkId][surfaceId] -> [u, v]

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

    void precalcApproxClosestSurfacePoint(int regionRadius, int precision = 0)
    {
        if (precision == 0) precision = precalcPrecision;
        approxClosestSurfacePoint.clear();

        for (int i = 0; i < surfaces.size(); i++)
        {
            const Surface & surface = surfaces[i];
            const auto validParams = surface.validParams(precision);

            for (const Chunk::ChunkId & chunkId : Chunk::ChunkId::chunkIdsInRadius(regionRadius))
            {
                Vector3D centre = chunkId.centre();
                long double minDist = DBL_MAX;
                std::pair<double, double> uv;
                for (auto [u, v] : validParams)
                {
                    if (long double dist = (surface(u, v) - centre).norm(); dist < minDist)
                    {
                        uv = std::make_pair(u, v);
                        minDist = dist;
                    }
                }
                approxClosestSurfacePoint[chunkId][i] = uv;
            }
        }
    }

    std::pair<double, double> preciseClosestSurfacePoint(const Vector3D & coord, const int & surfaceId, int precision = 0) const
    {
        if (precision == 0) precision = spreadPrecision;
        Chunk::ChunkId chunkId(coord);
        const Surface & surface = surfaces[surfaceId];
        
        std::pair<double, double> uv = approxClosestSurfacePoint.at(chunkId).at(surfaceId);
        long double minDist = (surface(uv.first, uv.second) - coord).norm();

        for (double duv = 0.125; duv > 0.01 / precision; )
        {
            auto [u0, v0] = uv;
            double u1 = u0 + double(1) / precision;
            double v1 = v0 + double(1) / precision;
            if (!surface.domain(u1, v0) || !surface.domain(u0, v1)) break;

            double gradU = minDist - (coord - surface(u1, v0)).norm();
            double gradV = minDist - (coord - surface(u0, v1)).norm();
            Vector3D /* Actually Vecotor2D */ grad = {gradU, gradV, 0};
            grad = grad.unit() * duv;
            double u2 = u0 + grad.x;
            double v2 = v0 + grad.y;
            if (!surface.domain(u2, v2))
            {
                duv /= 2;
                continue;
            }
            long double dist = (surface(u2, v2) - coord).norm();
            if (dist >= minDist)
            {
                duv /= 2;
                continue;
            }
            
            uv.first = u2, uv.second = v2;
            minDist = dist;
        }
        return uv;
    }

    Vector3D preciseNormalVectorAround(const Vector3D & coord, const int & surfaceId, int precision = 0) const
    {
        if (precision == 0) precision = spreadPrecision;
        Chunk::ChunkId chunkId(coord);
        const Surface & surface = surfaces[surfaceId];
        auto [u, v] = preciseClosestSurfacePoint(coord, surfaceId, precision);
        return surface.normalVectorAt(u, v, precision);
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

    std::list<std::pair<Vector3D, long double>> discreteSurfaceField(int precision, long double distance) const
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
                    intensity += freeCharge.q() * freeCharge.q() / deltaX.norm() / deltaX.norm() * deltaX.unit();
                }
                fields.emplace_back(extendedPoint, intensity.norm());
            }
        }
        return fields;
    }

    std::vector<std::vector<std::vector<std::optional<long double>>>>
        paramSurfaceField(int precision, long double distance, ElectricField outerfieldPlus = ElectricField())
        // paramSurfaceField[surfaceId][uInt][vInt] = intensity
    {
        const static long double K = 9E9;

        std::vector<std::vector<std::vector<std::optional<long double>>>> fields(surfaces.size());
        for (auto & row : fields)
        {
            row.resize(precision);
            for (auto & col : row) col.resize(precision);
        }

        #ifdef PERPENDICULAR_VERIFY
        int count = 0;
        long double cosOfAngle = 0.0;
        #endif

        for (int surfaceId = 0; surfaceId < surfaces.size(); surfaceId++)
        {
            const Surface & surface = surfaces[surfaceId];
            for (int uInt = 0; uInt < precision; uInt++) for (int vInt = 0; vInt < precision; vInt++)
            {
                double u = double(uInt) / precision, v = double(vInt) / precision;
                if (!surface.domain(u, v)) continue;
                
                Vector3D surfacePoint = surface(u, v);
                Vector3D normalVect = surface.normalVectorAt(u, v, precalcPrecision);
                if (normalVect == Vector3D::ZERO_VECTOR) continue;
                Vector3D extendedPoint = surfacePoint + distance * normalVect;
                Vector3D intensity;
                for (const FreeCharge & freeCharge : boundCharges)
                {
                    Vector3D deltaX = extendedPoint - freeCharge.coord;
                    intensity += K * freeCharge.q() / deltaX.norm() / deltaX.norm() * deltaX.unit();
                }

                intensity += outerfieldPlus.get(extendedPoint);

                // Gaussian curvature correction
                long double gaussR = sqrt(1 / surface.gaussianCurvatureAt(u, v, precalcPrecision));
                long double coef = (1 + (distance / gaussR)) * (1 + (distance / gaussR));

                fields[surfaceId][uInt][vInt] = intensity.norm() * coef;

                #ifdef PERPENDICULAR_VERIFY
                if (long double temp = intensity.cosineOfAngle(normalVect); !isnan(temp))
                {
                    count++;
                    cosOfAngle += temp;
                }
                #endif
            }
        }
        #ifdef PERPENDICULAR_VERIFY
        std::cout << "Averge cosOfAngle = " << cosOfAngle / count << '\n';
        #endif
        return fields;
    }

    std::vector<std::vector<std::vector<std::optional<long double>>>>
        paramSurfaceDensity(int precision, long double distance = 0.2 * Chunk::CHUNK_LENGTH)
        // paramSurfaceDensity[surfaceId][uInt][vInt] = count
    {
        std::vector<std::vector<std::vector<std::optional<long double>>>> density(surfaces.size());
        for (auto & row : density)
        {
            row.resize(precision);
            for (auto & col : row)
            {
                col.resize(precision);
                for (auto & element : col) element = 0.0;
            }
        }

        for (int surfaceId = 0; surfaceId < surfaces.size(); surfaceId++)
        {
            const Surface & surface = surfaces[surfaceId];
            for (const FreeCharge & charge : boundCharges)
            {
                const auto & [u, v] = preciseClosestSurfacePoint(charge.coord, surfaceId);
                Vector3D surfacePoint = surface(u, v);
                if ((surfacePoint - charge.coord).norm() >= distance) continue;
                int uInt = int(u * precision), vInt = int(v * precision);
                density[surfaceId][uInt][vInt] = *density[surfaceId][uInt][vInt] + 1;
            }
            for (int uInt = 0; uInt < precision; uInt++) for (int vInt = 0; vInt < precision; vInt++)
                density[surfaceId][uInt][vInt] = *density[surfaceId][uInt][vInt]
                    / surface.rawNormalVectorAt(double(uInt) / precision, double(vInt) / precision, precision).norm();
            // Divide with area of dudv
        }

        return density;
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
            double theta = u * PI, phi = v * 2 * PI;
            return centre + Vector3D{a * sin(theta) * cos(phi), b * sin(theta) * sin(phi), c * cos(theta)};
        };
        surface.domain = Surface::DOMAIN_FULL;
        surface.facingOut = true;

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
        return generateEllipse(centre, r, r, r, precalcPrecision, spreadPrecision);
    }

    static Conductor generateCube
    (
        const Vector3D & centre = Vector3D::ZERO_VECTOR, const double & l = 0.01, 
        int precalcPrecision = 128, int spreadPrecision = 64
    )
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

    static Conductor generateTimeglass
    (
        const Vector3D & centre = Vector3D::ZERO_VECTOR,
        const double & a = 0.04, const double & b = 0.04, const double & h = 0.06, 
        int precalcPrecision = 128, int spreadPrecision = 64
    )
    {
        Conductor conductor;

        Surface waist;
        waist.implicitSurface = [centre, a, b, h] (const Vector3D & coord) -> bool
        {
            const auto [x, y, z] = coord - centre;
            return (x / (a / 2)) * (x / (a / 2)) + (y / (b / 2)) * (y / (b / 2)) - (z / (h / 3) / sqrt(3)) * (z / (h / 3) / sqrt(3)) <= 1.0;
        };
        waist.paramSurface = [centre, a, b, h](const double & u, const double & v) -> Vector3D
        {
            double theta = v * 2 * PI, z = u * 2 * h - h;
            double c = sqrt(3) / 3 * h;
            return centre + Vector3D
            (
                a / 2 * sqrt(1 + z * z / c / c) * cos(theta), 
                b / 2 * sqrt(1 + z * z / c / c) * sin(theta),
                z
            );
        };
        waist.domain = Surface::DOMAIN_FULL;
        waist.facingOut = false;

        Surface upperBase = Surface::generatePlate(centre + Vector3D{-a, -b, h}, {2 * a, 0, 0}, {0, 2 * b, 0}, true, Surface::DOMAIN_CIRCLE);
        Surface lowerBase = Surface::generatePlate(centre + Vector3D{-a, -b, -h}, {2 * a, 0, 0}, {0, 2 * b, 0}, false, Surface::DOMAIN_CIRCLE);

        conductor.surfaces = { waist, upperBase, lowerBase };
        conductor.precalcPrecision = precalcPrecision;
        conductor.spreadPrecision = spreadPrecision;

        return conductor;
    }

    static Conductor generateTwistedTimeglass
    (
        const Vector3D & centre = Vector3D::ZERO_VECTOR,
        const double & a = 0.04, const double & b = 0.04, const double & h = 0.06, 
        int precalcPrecision = 128, int spreadPrecision = 64
    )
    {
        auto twist = [](const Vector3D & coord) -> Vector3D
        {
            const auto [x0, y0, z0] = coord;
            long double x = x0 + 128 * z0 * z0 * z0;
            long double y = y0 + 20 * x0 * x0 + 10 * z0 * z0;
            long double z = z0;
            return Vector3D{x, y, z};
        };

        auto untwist = [](const Vector3D & coord) -> Vector3D
        {
            const auto [x, y, z] = coord;
            long double z0 = z;
            long double x0 = x - 128 * z0 * z0 * z0;
            long double y0 = y - 20 * x0 * x0 - 10 * z0 * z0;
            
            return Vector3D{x0, y0, z0};
        };

        Conductor conductor;

        Surface waist;
        waist.implicitSurface = [centre, a, b, h, untwist] (const Vector3D & coord) -> bool
        {
            const auto [x, y, z] = untwist(coord - centre);
            return (x / (a / 2)) * (x / (a / 2)) + (y / (b / 2)) * (y / (b / 2)) - (z / (h / 3) / sqrt(3)) * (z / (h / 3) / sqrt(3)) <= 1.0;
        };
        waist.paramSurface = [centre, a, b, h, twist](const double & u, const double & v) -> Vector3D
        {
            double theta = v * 2 * PI, z0 = u * 2 * h - h;
            double c = sqrt(3) / 3 * h;
            long double x0 = a / 2 * sqrt(1 + z0 * z0 / c / c) * cos(theta);
            long double y0 = b / 2 * sqrt(1 + z0 * z0 / c / c) * sin(theta);
            return twist({x0, y0, z0}) + centre;
        };
        waist.domain = Surface::DOMAIN_FULL;
        waist.facingOut = false;

        auto generateTwistedPlate = [twist](const Vector3D & origin,
            const Vector3D & vecU, const Vector3D & vecV,
            bool facingOut = true, const Surface::Domain & domain = Surface::DOMAIN_FULL)
        {
            Surface plate;
            plate.paramSurface = [origin, vecU, vecV, twist]
                (const double & u, const double & v) -> Vector3D
            {
                return twist(origin + u * vecU + v * vecV);
            };
            plate.domain = domain;
            plate.facingOut = facingOut;
            plate.implicitSurface = [origin, vecU, vecV, facingOut]
                (const Vector3D & coord) -> bool
            {
                if (facingOut) return ((coord - origin) % (vecU * vecV) <= 0);
                return ((coord - origin) % (vecU * vecV) >= 0);
            };
            return plate;
        };

        Surface upperBase = generateTwistedPlate(centre + Vector3D{-a, -b, h}, {2 * a, 0, 0}, {0, 2 * b, 0}, true, Surface::DOMAIN_CIRCLE);
        Surface lowerBase = generateTwistedPlate(centre + Vector3D{-a, -b, -h}, {2 * a, 0, 0}, {0, 2 * b, 0}, false, Surface::DOMAIN_CIRCLE);

        conductor.surfaces = { waist, upperBase, lowerBase };
        conductor.precalcPrecision = precalcPrecision;
        conductor.spreadPrecision = spreadPrecision;

        return conductor;
    }

    static Vector3D xyzToRtp(const Vector3D & xyz) // r theta phi
    {
        auto & [x, y, z] = xyz;
        long double r = sqrt(x * x + y * y + z * z);
        long double theta = acos(z / r);
        long double phi = asin(y / (r * sin(theta)));
        return {r, theta, phi};
    }
    static Vector3D rtpToXyz(const Vector3D & rtp)
    {
        auto & [r, theta, phi] = rtp;
        return {r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta)};
    }

    static Conductor generateCthulhu
    (
        const Vector3D & centre = Vector3D::ZERO_VECTOR, 
        const double & size = 0.02, 
        int precalcPrecision = 128, int spreadPrecision = 128
    )
    {
        Conductor conductor;
        Surface surface;
        
        surface.implicitSurface = [centre, size](const Vector3D & coord) -> bool
        {
            const auto [r, theta, phi] = xyzToRtp(coord - centre);
            return r <= size * (2 + 0.1 * theta * theta * theta * sin(3 * phi) * sin(6 * theta));
        };
        surface.paramSurface = [size, centre](const double & u, const double & v) -> Vector3D
        {
            double theta = v * PI, phi = u * 2 * PI;
            double r = size * (2 + 0.1 * theta * theta * theta * sin(3 * phi) * sin(6 * theta));
            return centre + rtpToXyz({r, theta, phi});
        };
        surface.domain = Surface::DOMAIN_FULL;
        surface.facingOut = false;

        conductor.surfaces.push_back(surface);
        conductor.precalcPrecision = precalcPrecision;
        conductor.spreadPrecision = spreadPrecision;

        return conductor;
    }
};

const std::function<bool(const double &, const double &)>
    Conductor::Surface::DOMAIN_FULL = [](const double & u, const double & v)
    { return true; };
const std::function<bool(const double &, const double &)>
    Conductor::Surface::DOMAIN_CIRCLE = [](const double & u, const double & v)
    { return (u - 0.5) * (u - 0.5) + (v - 0.5) * (v - 0.5) <= 0.25; };

} // namespace HurryPeng

#endif
