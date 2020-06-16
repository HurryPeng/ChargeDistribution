#ifndef _HURRYPENG_CHUNK_HPP
#define _HURRYPENG_CHUNK_HPP

#include "Vector3D.hpp"
#include "ElectricField.hpp"
#include <initializer_list>
#include <set>
#include <array>
#include <string>
#include <sstream>

namespace HurryPeng
{

class Chunk
{
public:
    static const double CHUNK_LENGTH;

    struct ChunkId
    {
        int x;
        int y;
        int z;

        ChunkId() : x(0), y(0), z(0) {}
        ChunkId(const int & _x, const int &_y, const int &_z)
            :x(_x), y(_y), z(_z) {}
        ChunkId(const std::initializer_list<int> & il)
        {
            if (il.size() != 3) throw std::invalid_argument("Initializing a ChunkId with more ore less than 3 arguments");
            x = il.begin()[0], y = il.begin()[1], z = il.begin()[2];
        }

        ChunkId(const Vector3D & coord)
        {
            auto chunkFloor = [](const long double & x) -> int
            {
                return (int)(floor(x / CHUNK_LENGTH));
            };
            x = chunkFloor(coord.x);
            y = chunkFloor(coord.y);
            z = chunkFloor(coord.z);
        }

        bool operator==(const ChunkId & rhs) const
        {
            return x == rhs.x && y == rhs.y && z == rhs.z;
        }

        Vector3D centre() const
        {
            return CHUNK_LENGTH * Vector3D(x + 0.5, y + 0.5, z + 0.5);
        }

        std::array<Vector3D, 8> vertexes() const
        {
            std::array<Vector3D, 8> array;
            for (int i = 0; i < 8; i++)
                array[i] = CHUNK_LENGTH * Vector3D(x + (i & 1), y + ((i & 2) >> 1), z + ((i & 4) >> 2));
            return array;
        }

        bool covers(const Vector3D & coord) const
        {
           return *this == ChunkId(coord);
        }

        operator std::string() const
        {
            std::stringstream ss;
            ss << "{" << x << ", " << y << ", " << z << "}";
            return ss.str();
        }

        friend std::ostream & operator<<(std::ostream & lhs, const ChunkId & rhs);
    }; // struct ChunkId

    struct FreeCharge
    {
        const static long double Q; // Unit charge that all free chargrs share
        const static long double M; // Unit mass that all free charges share

        bool isPositive;
        Vector3D coord;
        Vector3D vel; // velocity

        FreeCharge(bool _isPositive, const Vector3D & _coord)
            :isPositive(_isPositive), coord(_coord) { }

        long double q() const { return isPositive ? Q : -Q; } // Returns signed q

        ChunkId chunkId() const { return ChunkId(coord); }
    }; // struct FreeCharge

private:
    ChunkId id;

    std::set<const FreeCharge *> freeCharges;
        // Pointers to FreeCharge objects subordinate to this Chunk
    
    bool fieldUpdated; // Whether chunkField is up-to-date
    ElectricField approxField;
    ElectricField preciseField;

public:
    Chunk() = delete;
    Chunk(const ChunkId & chunkId)
        :id(chunkId), fieldUpdated(true) {}

    ChunkId getId() const { return id; }
    Vector3D centre() const { return id.centre(); }
    std::array<Vector3D, 8> vertexes() const { return id.vertexes(); }
    bool covers(const Vector3D & coord) const { return id.covers(coord); }

    void assignCharge(const FreeCharge & charge)
    {
        #ifdef _DEBUG
        if (!covers(charge.coord))
            throw std::invalid_argument("Charge assigned to a Chunk not covering it");
        #endif
        freeCharges.insert(&charge);
        fieldUpdated = false;
    }
    void removeCharge(const FreeCharge & charge)
    {
        freeCharges.erase(&charge);
        fieldUpdated = false;
    }

    void updateStat()
    {
        if (fieldUpdated) return;

        long double positiveChargeSum = 0, negativeChargeSum = 0;
        int positiveChargeCount = 0, negativeChargeCount = 0;
        Vector3D positiveChargeCentre, negativeChargeCentre;
        for (const FreeCharge *charge : freeCharges)
        {
            if (charge->isPositive)
            {
                positiveChargeSum += charge->q();
                positiveChargeCount++;
                positiveChargeCentre += charge->coord;
            }
            else
            {
                negativeChargeSum += charge->q();
                negativeChargeCount++;
                negativeChargeCentre += charge->coord;
            }
            preciseField.overlayPointChargeField(charge->coord, charge->q());
        }
        approxField = ElectricField();
        if (positiveChargeCount != 0)
        {
            positiveChargeCentre /= positiveChargeCount;
            approxField.overlayPointChargeField(positiveChargeCentre, positiveChargeSum);
        }
        if (negativeChargeCount != 0)
        {
            negativeChargeCentre /= negativeChargeCount;
            approxField.overlayPointChargeField(negativeChargeCentre, negativeChargeSum);
        }
        fieldUpdated = true;
    }
    
    bool fieldIsUpdated() const { return fieldUpdated; }
    ElectricField getApproxField() const { return approxField; }
    ElectricField getPreciseField() const { return preciseField; }
};

const long double Chunk::FreeCharge::Q = 1.0E-12;
const long double Chunk::FreeCharge::M = 1.0;

const double Chunk::CHUNK_LENGTH = 0.01;

std::ostream & operator<<(std::ostream & lhs, const Chunk::ChunkId & rhs)
{
    lhs << std::string(rhs);
    return lhs;
}

} // namespace HurryPeng

#endif
