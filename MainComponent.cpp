#include "MeshRect.h"
#include <iomanip>
#include <fstream>

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

double u_0(double x, double y) {
    return exp(- x - y);
}

void getU_0(VectorXd& U_0) {
    double x, y;

    for (int n = 0; n < Nb; n++) {
        x = (n % (Nx + 1)) * h1;
        y = (n / (Nx + 1)) * h2;
        U_0(n) = u_0(x, y);
    }
}

int main(int argc, char *argv[]) {
    if (shape == "rectangle") {
        VectorXd U_0(Nb);

        getU_0(U_0);

        MeshRect mesh = MeshRect(U_0);
        mesh.applyBCtoU();

        for (int i = 0; i < MAX_ITERATIONS; i++) {
            cout << "ITERATION " << i+1 << endl;

            mesh.calculateAF();
            mesh.applyBCtoDelU();

            mesh.addDelU();

            if (mesh.endConditionMet()){
                break;
            }
        }

        mesh.output(argv[1]);
    }
}