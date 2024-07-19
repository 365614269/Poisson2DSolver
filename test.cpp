#include <iostream>
#include <fstream>
#include "globals.h"
#include "MeshRect.h"

using namespace std;

void output(string fileName) {
    ofstream fout("ans.vtk");

    fout << "# vtk DataFile Version 3.0" << endl;
    fout << "2D Poisson Equation Numeric Solution" << endl;
    fout << "ASCII" << endl;
    fout << "DATASET UNSTRUCTURED_GRID" << endl;
}

int main() {
    MeshRect mesh = MeshRect(U_0);
    long double x = 1.2;
    long double y = 0.8;

    cout << mesh.u(x, y) << endl;
}