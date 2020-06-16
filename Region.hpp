#ifndef _HURRYPENG_REGION_HPP
#define _HURRYPENG_REGION_HPP

#include "Vector3D.hpp"
#include "ElectricField.hpp"
#include "Chunk.hpp"
#include "Conductor.hpp"
#include <vector>
#include <list>
#include <map>
#include <initializer_list>
#include <algorithm>

namespace HurryPeng
{

class Region
{
private:
    const int radius = 8;
    std::map<Chunk::ChunkId, Chunk> chunks;
    std::vector<Conductor> conductors;
    ElectricField outerField;

    const long double dt;
    int tick;

public:
    Region() = delete;
    Region(const int & _radius, std::initializer_list<Conductor> _conductors,
        const ElectricField & _outerField, const long double & _dt)
        :radius(_radius), conductors(_conductors), outerField(_outerField), dt(_dt), tick(0)
    {
        for (const Chunk::ChunkId chunkId : Chunk::ChunkId::chunkIdsInRadius(radius))
            chunks.emplace(std::make_pair(chunkId, Chunk(chunkId)));
        
        for (Conductor & conductor : conductors)
        {
            conductor.precalcApproxNormalVector(radius);
            conductor.spreadCharges(true);

            for (const FreeCharge & charge : conductor.boundCharges)
                chunks.at(Chunk::ChunkId(charge.coord)).assignCharge(&charge);
        }

    }

    struct Stat
    {
        int tick = 0;
        int count = 0; // Number of free charges
        long double dxEstm = 0.0; // Average dx before elimination
        long double dxElim = 0.0; // Average dx after elimination
        long double vel = 0.0; // Average velocity
        long double acc = 0.0; // Average acceleration
        long double elimRate = 0.0; // Rate of eliminated charges
        long double bruteElimRate = 0.0; // Rate of cut-half-eliminated charges
        long double stayRate = 0.0; // Rate of charges that stayed
        long double superFoceRate = 0.0; // Rate of charges with excessive acceleration

        operator std::string() const
        {
            std::stringstream ss;
            ss << "Tick " << tick << " statistics: " << '\n';
            ss << "dxEstm, dxElim = " << dxEstm << ", " << dxElim << '\n';
            ss << "vel, acc = " << vel << ", " << acc << '\n';
            ss << "elim, brute = " << elimRate << ", " << bruteElimRate << '\n';
            ss << "stay, superforce = " << stayRate << ", " << superFoceRate << '\n';
            return ss.str();
        }
    };

    Stat proceed()
    {
        Stat stat;

        tick++;
        stat.tick = tick;
        for (auto & pr : chunks) pr.second.updateStat();

        for (Conductor & conductor : conductors)
        {
            for (FreeCharge & charge : conductor.boundCharges)
            {
                stat.count++;

                Vector3D acc = charge.accel(outerField);
                Chunk & curChunk = chunks.at(Chunk::ChunkId(charge.coord));
                for (const auto & pr : chunks)
                {
                    const Chunk & anyChunk = pr.second;
                    if (anyChunk.isSurrounding(curChunk))
                        acc += charge.accel(anyChunk.getPreciseField());
                    else acc += charge.accel(anyChunk.getApproxField());
                }

                // Max vel / acc limits
                if (charge.vel.norm() >= 0.01 / dt)
                    charge.vel = charge.vel.unit() * 0.01 / dt;
                if (acc.norm() >= 0.01 / dt / dt * 2)
                    acc = acc.unit() * 0.01 / dt / dt * 2;

                Vector3D dx = charge.vel * dt + acc * dt * dt / 2;
                if (dx.norm() >= 0.01) dx = dx.unit() * 0.01; // //

                if (acc.norm() >= 10000) stat.superFoceRate++, stat.acc -= acc.norm(); //std::cerr << acc << dx << '\n';
                stat.acc += acc.norm();
                charge.vel += acc * dt;

                stat.dxEstm += dx.norm();//

                std::vector<int> accumExceedIds = conductor.idOfSurfacesExceededBy(charge.coord + dx);
                if (accumExceedIds.size() == 1) // Exceeding a surface
                {
                    stat.elimRate++; //
                    const int & exceedId0 = accumExceedIds[0];
                    const Vector3D & norm0 = conductor.approxNormalVector[curChunk.getId()][exceedId0];

                    Vector3D dxPerpendicular = dx.projectOnto(norm0);
                    dx -= 1.0 * dxPerpendicular; // Erase perpendicular part from dx
                    charge.vel -= charge.vel.projectOnto(norm0);

                    if (!conductor.surfaces[accumExceedIds[0]](charge.coord + dx)) // Bounce
                        dx -= 0.5 * dxPerpendicular;

                    // Update accumeExceedIds to evaluate whether further operation is needed
                    for (const int & id : conductor.idOfSurfacesExceededBy(charge.coord + dx))
                        if (std::find(accumExceedIds.begin(), accumExceedIds.end(), id) == accumExceedIds.end())
                            accumExceedIds.push_back(id);
                }
                if (accumExceedIds.size() == 2) // Exceeding an edge
                {
                    const int & exceedId1 = accumExceedIds[0];
                    const int & exceedId2 = accumExceedIds[1];
                    const Vector3D & norm1 = conductor.approxNormalVector[curChunk.getId()][exceedId1];
                    const Vector3D & norm2 = conductor.approxNormalVector[curChunk.getId()][exceedId2];

                    Vector3D norm0 = (norm1 * norm2).unit();
                    dx = dx.projectOnto(norm0);
                    charge.vel = charge.vel.projectOnto(norm0);
                }
                // If no way has proved effective up till now, then we can only go to brute
                if (!conductor(charge.coord + dx))
                {
                    stat.bruteElimRate++;

                    for (int i = 1; i <= 8 && !conductor(charge.coord + dx); i++)
                        dx /= 2, charge.vel /= 2;
                    
                    if (!conductor(charge.coord + dx)) dx = charge.vel = Vector3D::ZERO_VECTOR;
                }

                if (dx == Vector3D::ZERO_VECTOR) stat.stayRate++;

                stat.vel += charge.vel.norm();
                stat.dxElim += dx.norm();

                charge.coord += dx;
                Chunk & nextChunk = chunks.at(Chunk::ChunkId(charge.coord));
                curChunk.removeCharge(&charge);
                nextChunk.assignCharge(&charge);
            }
        }

        stat.dxElim /= stat.count;
        stat.dxEstm /= stat.count;
        stat.vel /= stat.count;
        stat.acc /= stat.count;
        stat.elimRate /= stat.count;
        stat.bruteElimRate /= stat.count;
        stat.stayRate /= stat.count;
        stat.superFoceRate /= stat.count;

        return stat;
    }

    int getTick() const { return tick; }

    std::list<std::pair<Vector3D, long double>> calcSurfaceField(int precision, long double distance) const
    {
        std::list<std::pair<Vector3D, long double>> fields;
        for (const Conductor & conductor : conductors)
            fields.splice(fields.end(), conductor.calcSurfaceField(precision, distance));
        return fields;
    }

    std::list<FreeCharge> allCharges()
    {
        std::list<FreeCharge> freeCharges;
        for (const Conductor & conductor : conductors)
        {
            freeCharges.insert(freeCharges.end(),
                conductor.boundCharges.begin(), conductor.boundCharges.end());
        }
        return freeCharges;
    } 

}; // class Region

std::ostream & operator<<(std::ostream & lhs, const Region::Stat & rhs)
{
    lhs << std::string(rhs);
    return lhs;
}

} // namespace HurryPeng

#endif
