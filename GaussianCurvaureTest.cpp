#include "Region.hpp"
#include "DebugUtil.hpp"
#include "Timer.hpp"

using namespace std;
using namespace HurryPeng;

int main()
{
    Conductor sphere3 = Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.03, 128, 192);
    Conductor sphere6 = Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.06, 128, 192);
    Conductor ellipse = Conductor::generateEllipse(Vector3D::ZERO_VECTOR, 0.03, 0.03, 0.06, 128, 192);

    for (int uInt = 0; uInt < 64; uInt++)
    {
        cout << sphere3.surfaces[0].gaussianCurvatureAt((double)uInt / 64, 0, 128) << ' ';
    }
    cout << "\n\n";

    for (int uInt = 0; uInt < 64; uInt++)
    {
        cout << sphere6.surfaces[0].gaussianCurvatureAt((double)uInt / 64, 0, 128) << ' ';
    }
    cout << "\n\n";

    for (int uInt = 0; uInt < 64; uInt++)
    {
        cout << ellipse.surfaces[0].gaussianCurvatureAt((double)uInt / 64, 0, 128) << ' ';
    }
    cout << "\n\n";

    return 0;
}
