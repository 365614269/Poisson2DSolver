#include <iostream>
#include <utility>
#include "MeshRect.h"
#include "Node.h"
#include "integrate2D.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

long double lx = 4,ly = 3;
int Nx = 4, Ny = 3;

string shape = "rectangle";

long double f(long double x, long double y) {
    return 1;
}

long double delf(long double x, long double y) {
    return 0;
}

int main() {
    // Input for the attributes omitted

    if (shape == "rectangle") {
        MeshRect mesh = MeshRect(lx, ly, Nx, Ny, f, delf);
        mesh.calculateStiffness();
        mesh.displayStiffness();
    }
}