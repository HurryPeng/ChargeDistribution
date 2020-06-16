#define _DEBUG

#include <iostream>
#include <cmath>
#include "Chunk.hpp"
#include <list>
#include "DebugUtil.hpp"

using namespace std;
using namespace HurryPeng;

typedef Chunk::ChunkId ChunkId;

int main()
{
    //ChunkId id(Vector3D({0.655, 0.700, -0.301}));
    //for (const Vector3D & v : id.vertexes()) cout << v << '\n';
    //cout << id.centre() << '\n' << id.covers(Vector3D({0.655, 0.700, -0.301}));

    list<Chunk::FreeCharge> charges;
    charges.emplace_back(false, Vector3D{0.006, 0.007, 0});
    charges.emplace_back(true, Vector3D{0.007, 0.007, 0});
    charges.emplace_back(false, Vector3D{0.001, 0.003, 0});

    Chunk chunk({0, 0, 0});
    for (const Chunk::FreeCharge & charge : charges)
        chunk.assignCharge(&charge);
    chunk.updateStat();

    ElectricField field;
    field.overlayPointChargeField({0.006, 0.007, 0}, -Chunk::FreeCharge::Q);
    field.overlayPointChargeField({0.007, 0.007, 0}, +Chunk::FreeCharge::Q);
    field.overlayPointChargeField({0.001, 0.003, 0}, -Chunk::FreeCharge::Q);

    exportElectricField("out.txt", chunk.getPreciseField(), 10);

    return 0;
}
