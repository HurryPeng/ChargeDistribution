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
    const int radius = 8; // !!!!!!!!!!
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
            conductor.spreadCharges(true);///////////????

            if (conductor.boundCharges.empty()) std::cerr << "BoundCharges empy!\n"; //
            else std::cerr << "correct\n"; //

            for (const FreeCharge & charge : conductor.boundCharges)
                chunks.at(Chunk::ChunkId(charge.coord)).assignCharge(&charge);
        }

    }

    void proceed()
    {
        tick++;
        for (auto & pr : chunks) pr.second.updateStat();

        long double delta0 = 0, delta = 0;//
        long double avVel = 0, avAcc = 0;//
        //long double avAccApprox = 0;
        int eliminated = 0;//
        int superForce = 0;
        int brute = 0;
        int stayed = 0;
        int count = 0;//

        for (Conductor & conductor : conductors)
        {
            for (FreeCharge & charge : conductor.boundCharges)
            {
                count++;//
                //Vector3D accApprox = charge.accel(outerField);

                Vector3D acc = charge.accel(outerField);
                Chunk & curChunk = chunks.at(Chunk::ChunkId(charge.coord));
                for (const auto & pr : chunks)
                {
                    const Chunk & anyChunk = pr.second;
                    if (anyChunk.isSurrounding(curChunk))
                        acc += charge.accel(anyChunk.getPreciseField());
                    else acc += charge.accel(anyChunk.getApproxField());////////////////////
                    //if (anyChunk.isSurrounding(curChunk))//
                    //    accApprox += charge.accel(anyChunk.getPreciseField());//
                    //else accApprox += charge.accel(anyChunk.getApproxField());//
                }

                // Max vel / acc limits
                if (charge.vel.norm() >= 0.01 / dt) // //
                    charge.vel = charge.vel.unit() * 0.01 / dt;
                if (acc.norm() >= 0.01 / dt / dt * 2) // //
                    acc = acc.unit() * 0.01 / dt / dt * 2;

                Vector3D dx = charge.vel * dt + acc * dt * dt / 2;
                if (dx.norm() >= 0.01) dx = dx.unit() * 0.01; // //

                if (acc.norm() >= 10000) superForce++, avAcc -= acc.norm(); //std::cerr << acc << dx << '\n';
                avAcc += acc.norm();
                charge.vel += acc * dt;

                //Vector3D nextCoord = charge.coord + dx;
                delta0 += dx.norm();//

                std::vector<int> accumExceedIds = conductor.idOfSurfacesExceededBy(charge.coord + dx);
                //std::cerr << "Exceeds: " << accumExceedIds.size() << '\n';
                if (accumExceedIds.size() == 1) // Exceeding a surface
                {
                    eliminated++; //
                    const int & exceedId0 = accumExceedIds[0];
                    const Vector3D & norm0 = conductor.approxNormalVector[curChunk.getId()][exceedId0];

                    dx -= 1.25 * dx.projectOnto(norm0); // Erase perpendicular part from dx
                    // And only dx "bounces"! 
                    charge.vel -= charge.vel.projectOnto(norm0);

                    // Update accumeExceedIds to evaluate whether further operation is needed
                    for (const int & id : conductor.idOfSurfacesExceededBy(charge.coord + dx))
                        if (std::find(accumExceedIds.begin(), accumExceedIds.end(), id) == accumExceedIds.end())
                            accumExceedIds.push_back(id);
                }
                if (accumExceedIds.size() == 2) // Exceeding an edge
                {
                    //std::cerr << "Exceeding an edge\n";

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
                    brute++; //

                    // Prevent sticky superforce
                    if (dx.norm() > 10) dx = dx.unit() * 10;

                    for (int i = 1; i <= 8 && !conductor(charge.coord + dx); i++)
                        dx /= 2, charge.vel /= 2;
                    
                    if (!conductor(charge.coord + dx)) dx = charge.vel = Vector3D::ZERO_VECTOR;
                }
                
                //if (!conductor(charge.coord)) std::cerr << "???\n";
                /*

                if (!conductor(nextCoord))
                {
                    eliminated++;//
                    charge.vel = Vector3D::ZERO_VECTOR; //////////////////////
                    Vector3D & norm = conductor.approxNormalVector[curChunk.getId()];
                    nextCoord -= dx.projectOnto(norm); // Erase perpendicular part
                    for (int i = 1; i <= 4 && !conductor(nextCoord); i++)
                        nextCoord = (nextCoord + charge.coord) / 2;
                    if (!conductor(nextCoord)) nextCoord = charge.coord;
                }

                */

                if (dx == Vector3D::ZERO_VECTOR) stayed++;//

                avVel += charge.vel.norm();//
                //avAccApprox += accApprox.norm();//
                delta += dx.norm();//

                charge.coord += dx;
                try
                {
                    chunks.at(Chunk::ChunkId(charge.coord));
                }
                catch (std::exception & e)
                {
                    std::cerr << dx << '\n' << charge.coord << '\n';
                    std::abort();
                }
                Chunk & nextChunk = chunks.at(Chunk::ChunkId(charge.coord));
                curChunk.removeCharge(&charge);
                nextChunk.assignCharge(&charge);
            }
        }
        std::cerr << "Av. delta = " << delta0 / count << ", " << delta / count << '\n';
        std::cerr << "Av. vel, acc = " << avVel / count << ", " << avAcc / count << '\n';
        //std::cerr << "Av. accApprox = " << avAccApprox / count << '\n';
        std::cerr << "Eliminated, stayed = " << (double)eliminated / count  << ", " << (double)stayed / count << '\n';
        std::cerr << "Superforce, brute = " << (double)superForce / count << ", " << (double)brute / count << '\n';
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

} // namespace HurryPeng

#endif
