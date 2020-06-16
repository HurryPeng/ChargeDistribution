#include "ClosedSurface.hpp"
#include "Vector3D.hpp"
#include "DebugUtil.hpp"
#include <cmath>

using namespace std;
using namespace HurryPeng;

int main()
{
    exportPointSet("points.txt", ClosedSurface::generateCube({1, 0, 2}, 1).surfacePoints(100));

    return 0;
}
