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
    bool initialised = false;

    const int radius = 8;
    std::map<Chunk::ChunkId, Chunk> chunks;
    std::vector<Conductor> conductors;

    long double dt;
    long double elapsed = 0.0;

    int tick;
    long double dynamicConvergenceTime;
        // If set to 0.0, then dynamic tick adjustment will be disabled. 
        // Otherwise, dt will be adjusted every tick to make average dx
        // begin with 0.8 * CHUNK_LENGTH on tick 0, and gradually decay to
        // 0.2 * CHUNK_LENGTH when elpased >= dynamicConvergenceTick

public:
    ElectricField outerField;

    std::string summary;

    Region() = delete;
    Region(const int & _radius, std::initializer_list<Conductor> _conductors, 
        const ElectricField & _outerField, const long double & _dt, 
        const std::string & _summary = "Untitled", const long double & _dynamicConvergenceTime = 0.0, 
        const std::list<Vector3D> & _importedPointSet = std::list<Vector3D>())
        :radius(_radius), conductors(_conductors), outerField(_outerField), dt(_dt),
        tick(0), summary(_summary), dynamicConvergenceTime(_dynamicConvergenceTime)
        {
            if (!_importedPointSet.empty()) initialise(_importedPointSet);
        }

    void initialise(const std::list<Vector3D> & importedPointSet = std::list<Vector3D>())
    {
        initialised = true;

        for (const Chunk::ChunkId chunkId : Chunk::ChunkId::chunkIdsInRadius(radius))
            chunks.emplace(std::make_pair(chunkId, Chunk(chunkId)));
        
        for (Conductor & conductor : conductors)
        {
            conductor.precalcApproxClosestSurfacePoint(radius);
            if (importedPointSet.empty()) conductor.spreadCharges(true);
            else
            {
                for (const Vector3D & point : importedPointSet)
                    if (conductor(point)) conductor.boundCharges.emplace_back(true, point);
            }

            for (const FreeCharge & charge : conductor.boundCharges)
                chunks.at(Chunk::ChunkId(charge.coord)).assignCharge(&charge);
        }
    }

    struct Stat
    {
        int tick = 0;
        int count = 0; // Number of free charges
        long double dt = 0.0; // Current dt
        long double elapsed = 0.0; // Elapsed time
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
            ss << "dt, elapsed = " << dt << ", " << elapsed << '\n';
            ss << "dxEstm, dxElim = " << dxEstm << ", " << dxElim << '\n';
            ss << "vel, acc = " << vel << ", " << acc << '\n';
            ss << "elim, brute = " << elimRate << ", " << bruteElimRate << '\n';
            ss << "stay, superforce = " << stayRate << ", " << superFoceRate << '\n';
            return ss.str();
        }
    };

    Stat proceed()
    {
        if (!initialised) initialise();

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
                for (const auto & [chunkId, chunk] : chunks)
                {
                    if (chunk.isSurrounding(curChunk))
                        acc += charge.accel(chunk.getPreciseField());
                    else acc += charge.accel(chunk.getApproxField());
                }

                // Max vel / acc limits
                if (charge.vel.norm() >= Chunk::CHUNK_LENGTH / dt)
                    charge.vel = charge.vel.unit() * Chunk::CHUNK_LENGTH / dt;
                if (acc.norm() >= Chunk::CHUNK_LENGTH / dt / dt * 2)
                    acc = acc.unit() * Chunk::CHUNK_LENGTH / dt / dt * 2;

                // velocity is "nerfed" to 0.7 to avoid too much oscillation
                charge.vel *= 0.7;
                Vector3D dx = charge.vel * dt + acc * dt * dt / 2;
                if (dx.norm() >= Chunk::CHUNK_LENGTH) dx = dx.unit() * Chunk::CHUNK_LENGTH; // //

                if (acc.norm() >= 1E8 * Chunk::CHUNK_LENGTH * Chunk::CHUNK_LENGTH)
                    stat.superFoceRate++, stat.acc -= acc.norm(); //std::cerr << acc << dx << '\n';
                stat.acc += acc.norm();
                charge.vel += acc * dt;

                stat.dxEstm += dx.norm();//

                std::vector<int> accumExceedIds = conductor.idOfSurfacesExceededBy(charge.coord + dx);
                if (accumExceedIds.size() == 1) // Exceeding a surface
                {
                    stat.elimRate++; //
                    const int & exceedId0 = accumExceedIds[0];
                    const Vector3D & norm0 = conductor.preciseNormalVectorAround(charge.coord, exceedId0);
                    
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
                    const Vector3D & norm1 = conductor.preciseNormalVectorAround(charge.coord, exceedId1);
                    const Vector3D & norm2 = conductor.preciseNormalVectorAround(charge.coord, exceedId2);

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

        elapsed += dt;

        stat.dt = dt;
        stat.elapsed = elapsed;
        stat.dxElim /= stat.count;
        stat.dxEstm /= stat.count;
        stat.vel /= stat.count;
        stat.acc /= stat.count;
        stat.elimRate /= stat.count;
        stat.bruteElimRate /= stat.count;
        stat.stayRate /= stat.count;
        stat.superFoceRate /= stat.count;

        if (dynamicConvergenceTime != 0.0)
        {
            long double dxExpected = 0.2 * Chunk::CHUNK_LENGTH;
            if (elapsed <= dynamicConvergenceTime)
                dxExpected = (0.5 + 0.3 * cos(PI / dynamicConvergenceTime * elapsed)) * Chunk::CHUNK_LENGTH;
            dt = dxExpected / stat.dxEstm * dt;
        }

        return stat;
    }

    int getTick() const { return tick; }

    const std::vector<Conductor> & getConductors() const { return conductors; }

    std::list<std::pair<Vector3D, long double>> discreteSurfaceField(int precision, long double distance) const
    {
        std::list<std::pair<Vector3D, long double>> fields;
        for (const Conductor & conductor : conductors)
            fields.splice(fields.end(), conductor.discreteSurfaceField(precision, distance));
        return fields;
    }

    std::vector<std::vector<std::vector<std::vector<std::optional<long double>>>>>
        paramSurfaceField(int precision, long double distance)
        // paramSurfaceField[conductorId][surfaceId][uInt][vInt] = intensity
    {
        std::vector<std::vector<std::vector<std::vector<std::optional<long double>>>>> fields(conductors.size());
        for (int conductorId = 0; conductorId < conductors.size(); conductorId++)
            fields[conductorId] = conductors[conductorId].paramSurfaceField(precision, distance, outerField);

        return fields;
    }

    std::vector<std::vector<std::vector<std::vector<std::optional<long double>>>>>
        paramSurfaceDensity(int precision, long double distance = 0.2 * Chunk::CHUNK_LENGTH)
        // paramSurfaceDensity[conductorId][surfaceId][uInt][vInt] = count
    {
        std::vector<std::vector<std::vector<std::vector<std::optional<long double>>>>> density(conductors.size());
        for (int conductorId = 0; conductorId < conductors.size(); conductorId++)
            density[conductorId] = conductors[conductorId].paramSurfaceDensity(precision, distance);

        return density;
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
