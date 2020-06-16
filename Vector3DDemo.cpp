#include <iostream>
#include "Vector3D.hpp"

using namespace std;
using namespace HurryPeng;

int main()
{
    Vector3D v1 = {2, 0, 0}, v2 = {1.6, 2, 0};
    cout << v1 << ' ' << v2 << ' ' << v1.outerProduct(v2) << ' ' << v1.innerProduct(v2) << '\n';

    Vector3D v3 = v1 + v2;
    cout << v3.x << ' ' << v3.unit() << ' ' << v3.norm() << '\n';

    return 0;
}
