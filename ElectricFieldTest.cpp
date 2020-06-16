#include <iostream>
#include "Vector3D.hpp"
#include "ElectricField.hpp"
#include <fstream>
#include "DebugUtil.hpp"

using namespace std;
using namespace HurryPeng;

int main()
{
    ofstream ofs("out.txt");

    ElectricField electricField;

    electricField.overlayUniformField({1, 1, 1});
    electricField.overlayPointChargeField({-1, 0, 0}, 5E-9);
    electricField.overlayPointChargeField({1, 0, 0}, -10E-9);
    // electricField.overlayPointChargeField({2, 3, 1}, 10E-9);
    electricField.overlayPointChargeField({-1, -3 ,0}, -5E-9);

    ofs << "{\n";
    for (int i = -50; i < 50; i++) for (int j = -50; j < 50; j++)
    {
        Vector3D temp = electricField.get({double(i) / 10, double(j) / 10, 0});
        ofs << "    {{" << double(i) / 10 << ", " << double(j) / 10 << "}, {" << temp.x << ", " << temp.y << "}}";
        if (i != 49 || j != 49) ofs << ',';
        ofs << '\n';
    }
    ofs << "}\n";

    return 0;
}
