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
    //Conductor conductor = Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.06, 64, 64);
    Conductor conductor = Conductor::generateEllipse(Vector3D::ZERO_VECTOR, 0.03, 0.03, 0.06, 128, 128);
    ElectricField outerField;
    //outerField.overlayUniformField(Vector3D{2E4, 0, 0});

    //conductor.spreadCharges(true);///////////////
    //exportPointSet("points.txt", chargesToPoints(conductor.boundCharges));

    const double TICK_SECOND = 0.0008;
    Region region(8, {conductor}, outerField, TICK_SECOND);

    cout << "Precalc finished " << t1.elapsed() << '\n';

    if (region.allCharges().empty()) throw "FUCK";

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
        for (const pair<Vector3D, long double> & pr : region.calcSurfaceField(8, 0.006))
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
            region.proceed();

            cout << "Tick " << region.getTick() << ", " << region.getTick() * TICK_SECOND << " seconds finished " << t1.elapsed() << "\n\n";
        }
        cout << "group finished \n";
        //getchar();
    }

    return 0;
}
