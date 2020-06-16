#include "Conductor.hpp"
#include "Vector3D.hpp"
#include "DebugUtil.hpp"
#include <cmath>

using namespace std;
using namespace HurryPeng;

int main()
{
    auto conductor = Conductor::generateCube({0, 0, 0}, 0.06, 128, 64);
    conductor.precalcApproxClosestSurfacePoint(8);
    conductor.spreadCharges(true);

    //for (auto & v : temp) cout << v.norm() << ' ';
    exportPointSet("points.txt", chargesToPointSet(conductor.boundCharges));

    cout << "1";
    /*
    ofstream ofs("normvec.txt");
    conductor.precalcApproxNormalVector(8, 100);
    //const std::pair<const HurryPeng::Chunk::ChunkId, HurryPeng::Vector3D> * prev = nullptr;
    cout << "2";
    for (auto & pr : conductor.approxNormalVector)
    {
        ofs << "{" << pr.first.centre() << ", " << pr.second << "}, \n";
    }
    ofs.close();
    */

    return 0;
}
