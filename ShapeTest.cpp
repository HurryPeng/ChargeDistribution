#include "Conductor.hpp"
#include "Vector3D.hpp"
#include "DebugUtil.hpp"
#include <cmath>

using namespace std;
using namespace HurryPeng;

int main()
{
    Conductor conductor = Conductor::generateTimeglass(Vector3D::ZERO_VECTOR, 0.04, 0.04, 0.06, 128, 64);
    conductor.spreadCharges(true);

    exportPointSet("points.txt", chargesToPointSet(conductor.boundCharges));

    return 0;
}
