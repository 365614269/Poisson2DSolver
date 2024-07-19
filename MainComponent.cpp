#include "MeshRect.h"
#include <iomanip>
#include <fstream>

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

bool EC(VectorXd delU, VectorXd U, VectorXd F) {
    long double abserr = delU.norm() / U.norm();
    long double relerr = F.norm();
    cout << "Abs. Error: " << abserr << " ;;; Rel. Error: " << relerr << endl;

    return (abserr < ABSOLUTE_TOLERANCE) && (relerr < RELATIVE_TOLERANCE);
}

int main() {
    if (shape == "rectangle") {
        ofstream fout("ans.txt");
        MeshRect mesh = MeshRect(U_0);
        mesh.applyBCtoU();

        for (int i = 0; i < MAX_ITERATIONS; i++) {
            cout << "ITERATION " << i+1 << endl;

            mesh.calculateAF();
            mesh.applyBCtoDelU();

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

        fout << U << endl;
    }
}