#include <iostream>
#include <utility>
#include <cmath>
#include "MeshRect.h"
#include "Node.h"
#include "integrate2D.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

long double lx = 4,ly = 3;
int Nx = 4, Ny = 3;

#define MAX_ITERATIONS 10

string shape = "rectangle";


long double f(long double x, long double y) {
    return 1;
}

long double delf(long double x, long double y) {
    return 0;
}

bool EC(VectorXd delU, VectorXd U, VectorXd F) {
    return (delU.norm() / U.norm() < 10e-4) && (F.norm() < 10e-4);
}

int main() {
    // Input for the attributes omitted

    if (shape == "rectangle") {
        VectorXd U_0 = VectorXd::Ones((Nx + 1) * (Ny + 1));
        MeshRect mesh = MeshRect(lx, ly, Nx, Ny, f, delf, U_0);

        for (int i = 0; i < MAX_ITERATIONS; i++) {
            mesh.calculateAF();
            mesh.applyBC();
            MatrixXd A = mesh.getStiffness();
            VectorXd F = mesh.getF();
            VectorXd delU = A.colPivHouseholderQr().solve(F);
            VectorXd U = mesh.getU();
            VectorXd U_new = U + delU;
            mesh.setU(U_new);

            if (EC(delU, U, F)) {
                break;
            }
        }

        VectorXd U = mesh.getU();

        cout << U << endl;
    }
}