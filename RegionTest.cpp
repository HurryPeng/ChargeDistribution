#define _DEBUG

#include "Region.hpp"
#include "DebugUtil.hpp"
#include <boost/timer.hpp>

using namespace std;
using namespace HurryPeng;
using namespace boost;

int main()
{
    timer t1;
    Conductor conductor = Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.06, 128, 192); // Case 1, 6
    //Conductor conductor = Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.06, 128, 64); // Case 0
    //Conductor conductor = Conductor::generateCube(Vector3D::ZERO_VECTOR, 0.12, 64, 96); // Case 5
    //Conductor conductor = Conductor::generateEllipse(Vector3D::ZERO_VECTOR, 0.03, 0.03, 0.06, 128, 192); // Case 2
    //Conductor conductor = Conductor::generateTimeglass(Vector3D::ZERO_VECTOR, 0.04, 0.04, 0.06, 128, 64); // Case 3
    //Conductor conductor = Conductor::generateTimeglass(Vector3D::ZERO_VECTOR, 0.04, 0.04, 0.06, 128, 128); // Case 4

    ElectricField outerField;
    //outerField.overlayUniformField(Vector3D{0, 0, 4E4}); // Case 1
    //outerField.overlayUniformField(Vector3D{0, 0, 4E3}); // Case 0
    //outerField.overlayPointChargeField(Vector3D{0.04, 0.04, 0.06}, 2.4E-8); // Case 4
    //outerField.overlayPointChargeField(Vector3D{0.03, 0, 0}, -1.6E-8); // Case 4
    outerField.overlayPointChargeField(Vector3D{0, 0, 0.12}, -1.0E-7); // Case 6

    //conductor.spreadCharges(true);///////////////
    //exportPointSet("points.txt", chargesToPoints(conductor.boundCharges));

    //const double TICK_SECOND = 0.004; // Case 5
    //const double TICK_SECOND = 0.003; // Case 0
    const double TICK_SECOND = 0.002; // Case 6
    //const double TICK_SECOND = 0.001; // Case 1, 2, 4, 6
    Region region(8, {conductor}, outerField, TICK_SECOND);

    cout << "Precalc finished " << t1.elapsed() << '\n';

    while (true)
    {
        list<Vector3D> toPrint = chargesToPoints(region.allCharges());
        for (const auto & charge : region.allCharges())
        {
            if (charge.vel.norm() >= 0.01 / 0.002) toPrint.push_back(charge.coord);
        }
        ofstream ofs("result.txt");
        double min = DBL_MAX, max = DBL_MIN;
        Vector3D vctMin, vctMax;
        for (const pair<Vector3D, long double> & pr : region.calcSurfaceField(64, 0.003))
        {
            if (pr.second < min) min = pr.second, vctMin = pr.first;
            if (pr.second > max) max = pr.second, vctMax = pr.first;
            ofs << pr.second << '\n';
            //toPrint.push_back(pr.first);
        }
        //toPrint.push_back(vctMin);
        //toPrint.push_back(vctMax);
        exportPointSet("points.txt", toPrint);
        cout << "min max phi = " << min << ", " << max << "\n";
        cout << vctMin << " and " << vctMax << "\n\n";

        for (int i = 1; i <= 1; i++)
        {
            t1.restart();
            cout << region.proceed();

            cout << "Tick " << region.getTick() << ", " << region.getTick() * TICK_SECOND << " seconds finished " << t1.elapsed() << "\n\n";
        }
        cout << "group finished \n";
        getchar();
    }

    return 0;
}
