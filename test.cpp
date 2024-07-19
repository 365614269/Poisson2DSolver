#include <iostream>
#include "globals.h"
#include "MeshRect.h"

using namespace std;

int main() {
    MeshRect mesh = MeshRect(U_0);
    long double x = 1.2;
    long double y = 0.8;

    cout << mesh.u(x, y) << endl;
}