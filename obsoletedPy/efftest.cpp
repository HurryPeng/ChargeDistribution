#include <cstdio>
#include <cmath>

using namespace std;

// long long long a;

void loop(int &a)
{
    for (int j = 1; j <= 1000000; j++) a += sqrt(2);
}

int main()
{
    int a = 0;
    for (int i = 1; i <= 1000; i++)
    {
        loop(a);
        //printf("%d\n", i);
    }
    printf("%d finished\n", a);
    return 0;
}
