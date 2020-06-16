// ChargeDistributionDemo.cpp
// HurryPeng

#include "Region.hpp"
#include "DebugUtil.hpp"
#include "Timer.hpp"

using namespace std;
using namespace HurryPeng;

// Preset testcases

Region testcase0
(
    8,
    { Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.06, 128, 64) },
    ElectricField(),
    0.003,
    "Sphere Conductor, Precision 64"
);

Region testcase1
(
    8,
    { Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.06, 128, 192) },
    ElectricField(),
    0.001,
    "Sphere Conductor, Precision 192"
);

Region testcase2
(
    8,
    { Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.06, 128, 192) },
    ElectricField().overlayUniformField({0, 0, 2E4}),
    0.001,
    "Sphere Conductor in Uniform Field, Precision 192", 
    0.16
);

Region testcase3
(
    8,
    { Conductor::generateSphere(Vector3D::ZERO_VECTOR, 0.06, 128, 192) },
    ElectricField().overlayPointChargeField(Vector3D{0, 0, 0.24}, -2E-7),
    0.001,
    "Sphere Conductor in Point Charge Field, Precision 192", 
    0.16
);

Region testcase4
(
    8,
    { Conductor::generateEllipse(Vector3D::ZERO_VECTOR, 0.03, 0.03, 0.06, 128, 192) },
    ElectricField(),
    0.001,
    "Ellipse Conductor, Precision 192", 
    0.16
);

Region testcase5
(
    8,
    { Conductor::generateCube(Vector3D::ZERO_VECTOR, 0.12, 64, 96) },
    ElectricField(),
    0.004,
    "Cube Conductor, Precision 96", 
    0.16
);

Region testcase6
(
    8,
    { Conductor::generateTimeglass(Vector3D::ZERO_VECTOR, 0.04, 0.04, 0.06, 128, 64) },
    ElectricField(),
    0.004,
    "Timeglass Conductor, Precision 64", 
    0.16
);

Region testcase7
(
    8,
    { Conductor::generateTimeglass(Vector3D::ZERO_VECTOR, 0.04, 0.04, 0.06, 128, 128) },
    ElectricField(),
    0.001,
    "Timeglass Conductor, Precision 128", 
    0.16
);

Region testcase8
(
    8,
    { Conductor::generateTimeglass(Vector3D::ZERO_VECTOR, 0.04, 0.04, 0.06, 128, 128) },
    ElectricField()
        .overlayPointChargeField(Vector3D{0.04, 0.04, 0.06}, 2.4E-8)
        .overlayPointChargeField(Vector3D{0.03, 0, 0}, -1.6E-8), 
    0.001,
    "Timeglass Conductor in Two Point Charge Fields, Precision 128", 
    0.16
);

Region testcase9
(
    8,
    { Conductor::generateTwistedTimeglass(Vector3D{0, -0.03, 0}, 0.04, 0.04, 0.06, 128, 64) },
    ElectricField(), 
    0.004,
    "Twisted Timeglass Conductor, Precision 64", 
    0.16
);

Region testcase10
(
    8,
    { Conductor::generateTwistedTimeglass(Vector3D{0, -0.03, 0}, 0.04, 0.04, 0.06, 128, 128) },
    ElectricField(), 
    0.001,
    "Twisted Timeglass Conductor, Precision 128", 
    0.16
);

Region testcase11
(
    8,
    { Conductor::generateTwistedTimeglass(Vector3D{0, -0.03, 0}, 0.04, 0.04, 0.06, 128, 128) },
    ElectricField()
        .overlayPointChargeField(Vector3D{0.04, 0.06, 0.06}, 2.4E-8)
        .overlayPointChargeField(Vector3D::ZERO_VECTOR, -1.6E-8),
    0.001,
    "Twisted Timeglass Conductor in Two Point Charge Fields, Precision 128", 
    0.16
);

Region testcase12
(
    8,
    { Conductor::generateCthulhu(Vector3D{0, 0, 0.02}, 0.02, 128, 192) },
    ElectricField(),
    0.001,
    "Cthulhu Conductor, Precision 192", 
    0.16
);

Region testcase13
(
    8,
    {
        Conductor::generateSphere({0, 0, 0.03}, 0.03, 128, 168),
        Conductor::generateSphere({0, 0, -0.03}, 0.03, 128, 168),
    },
    ElectricField(),
    0.001,
    "Two Sphere Conductors, Precision 168", 0.16
);

int main()
{
    Region & region = testcase4; // Select testcase here
    cout << region.summary << '\n';
    Timer timer;
    region.initialise();
    cout << "Precalculation finished in " << timer.read() << " seconds. \n\n";

    while (region.getTick() <= 128) // Emulation will keep running until manually stopped
    {
        cout << "Exporting point set... ";
        exportPointSet("points.txt", chargesToPointSet(region.allCharges()));
            // points.txt is adapted for Mathematica plotting
            // Run the following Wolfram languange code in Mathematica (or anything Wolfram) to draw a plot: 
            // pointfile = FileNameJoin[{NotebookDirectory[], "points.txt"}];
            // ListPointPlot3D[ToExpression[Import[pointfile]], Axes -> True, BoxRatios -> Automatic]
        // exportPointSetForPy("pointsPy.txt", chargesToPointSet(region.allCharges()));
            // If you don't have Mathematica, Python with matplotlib can be a good alternative
            // Run plot.py, and then it will search for pointsPy.txt and plot it
            // However, it will be way slower, which is not surprising at all when using anything Python
        cout << "done\n\n";

        bool doStat = true;
        if (doStat)
        {
            cout << "Distribution statistics: \n";
            for (const long double & intensity : paramStatU(region.paramSurfaceField(64, 0.001)[0][0]))
                cout << intensity << ' ';
            cout << "\n\n";
        }

        int group = 1; // Export points and do statistics every `group` ticks
        for (int i = 0; i < group; i++)
        {
            timer.restart();
            cout << region.proceed();
            cout << "Tick " << region.getTick() << " finished in " << timer.read() << " seconds. \n\n";
        }
        cout << "Group finished\n";
    }

    cin.get();
    cin.get();

    return 0;
}
