#include "MeshRect.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

MeshRect::MeshRect(VectorXd U_) {
    this->Uv = U_;
    this->Fv = VectorXd::Zero(Nb);
    this->stiffness = MatrixXd::Zero(Nb, Nb);

    for (int i = 0; i < Ne; i++) {
        this->elements[i] = Element(i);
    }
}

Element MeshRect::getElement(int n) {
    return this->elements[n];
}

int MeshRect::Tb(int localNodeIndex, int i, int j) {
    int offset[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    int x = i + offset[localNodeIndex][0];
    int y = j + offset[localNodeIndex][1];

    return x * (Nx + 1) + y;  // exchange index
}

int MeshRect::Tb(int localNodeIndex, int n) {
    int offset[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    pair<int, int> pos = make_pair(n / Nx, n % Nx);  // exchange index

    int x = pos.first + offset[localNodeIndex][0];
    int y = pos.second + offset[localNodeIndex][1];

    return x * (Nx + 1) + y;  // exchange index
}

void MeshRect::addAF(int elemIndex) {
    Element element = this->getElement(elemIndex);
    long double ans1, ans2, ans3;
    long double x1,x2,y1,y2;

    pair<int, int> index = make_pair(elemIndex / Nx, elemIndex % Nx);  // exchange index
    x1 = index.second * h1;
    x2 = (index.second + 1) * h1;
    y1 = (Ny - index.first - 1) * h2;
    y2 = (Ny - index.first) * h2;

    for (int alpha = 0; alpha < Nlb; alpha++) {
        ans3 = 0;

        for (int beta = 0; beta < Nlb; beta++) {
            ans1 = this->integrateDelpsi(element.nodes[beta], element.nodes[alpha], x1, x2, y1, y2);
            ans2 = this->integratePsiDelF(element.nodes[alpha], element.nodes[beta], x1, x2, y1, y2);
            ans3 += ans1 * this->Uv(Tb(beta, elemIndex));

            long double ans = -ans1-ans2;
            this->stiffness(Tb(beta, elemIndex), Tb(alpha, elemIndex)) += ans;
        }

        ans3 += this->integrateFPsi(element.nodes[alpha], x1, x2, y1, y2);
        this->Fv(Tb(alpha, elemIndex)) += ans3;
    }
}

void MeshRect::calculateAF() {
    this->stiffness = MatrixXd::Zero(Nb, Nb);
    this->Fv = VectorXd::Zero(Nb);

    for (int i = 0; i < Ne; i++) {
        this->addAF(i);
    }
}

MatrixXd MeshRect::getStiffness() {
    return this->stiffness;
}

VectorXd MeshRect::getF() {
    return this->Fv;
}

VectorXd MeshRect::getU() {
    return this->Uv;
}

void MeshRect::setU(VectorXd U_n) {
    this->Uv = U_n;
}

void MeshRect::applyBCtoU() {
    for (int n = 0; n < Nb; n++) {
        pair<int, int> index = make_pair(n / (Nx + 1), n % (Nx + 1));  // exchange index

        if (index.second == 0 || index.second == Nx) {
            long double x = index.second * h1;
            long double y = Ny - index.first * h2;

            if (this->g(x, y) != 999.9) {
                this->Uv(n) = this->g(x, y);
            }
        }
    }
}

void MeshRect::applyBCtoDelU() {
    for (int n = 0; n < Nb; n++) {
        pair<int, int> index = make_pair(n / (Nx + 1), n % (Nx + 1));  // exchange index

        if (index.second == 0 || index.second == Nx) {
            for (int i = 0; i < Nb; i++) {
                this->stiffness(n,i) = 0;
            }
            this->stiffness(n,n) = 1;
            this->Fv(n) = 0;
        }
    }
}

pair<int, int> MeshRect::elem(long double x, long double y) {
    int a = (ly - y) / h2;
    int b = x / h1;

    a = min(Ny, a);
    b = min(Nx, b);

    return make_pair(a, b);
}

long double MeshRect::u(long double x, long double y){
    //find which element (x, y)is in.
    pair<int, int> index = this->elem(x, y);
    Element element = this->getElement(index.first * Nx + index.second);  // exchange index

    VectorXd v1(Nlb);
    VectorXd v2(Nlb);

    // add up the 4 components * shape functions of the 4 nodes

    for (int i = 0; i < Nlb; i++) {
        v1(i) = element.nodes[i].psi(x, y);
        int n = Tb(i, index.first, index.second);
        v2(i) = this->Uv(n);
    }

    return v1.dot(v2);
}

long double MeshRect::f(long double x, long double y) {
    return exp(-this->u(x, y));
    // return 0;
}

long double MeshRect::delf(long double x, long double y) {
    return -exp(-this->u(x, y));
    // return 0;
}

long double MeshRect::g(long double x, long double y) {
    if (x == 0) {
        return 1;
    } else if (x == lx) {
        return 2;
    } else {
        return 999.9;
    }
}

long double MeshRect::integratePsiDelF(Node node1, Node node2, long double a1, long double b1, long double a2, long double b2) {
    long double diff1 = (b1 - a1) / 2;
    long double avg1 = (b1 + a1) / 2;
    long double diff2 = (b2 - a2) / 2;
    long double avg2 = (b2 + a2) / 2;

    long double ans = 0;
    long double x,y;

    for (int i = 0; i < PRECISION; i++) {
        for (int j = 0; j < PRECISION; j++) {
            x = diff1 * abscissae[i] + avg1;
            y = diff2 * abscissae[j] + avg2;
            ans += weights[i] * weights[j] * node1.psi(x, y) * node2.psi(x, y) * this->delf(x, y);
        }
    }

    return diff1 * diff2 * ans;
}

long double MeshRect::integrateDelpsi(Node node1, Node node2, long double a1, long double b1, long double a2, long double b2) {
    long double diff1 = (b1 - a1) / 2;
    long double avg1 = (b1 + a1) / 2;
    long double diff2 = (b2 - a2) / 2;
    long double avg2 = (b2 + a2) / 2;

    long double ans = 0;
    long double x,y;

    for (int i = 0; i < PRECISION; i++) {
        for (int j = 0; j < PRECISION; j++) {
            x = diff1 * abscissae[i] + avg1;
            y = diff2 * abscissae[j] + avg2;
            ans += weights[i] * weights[j] * node1.delpsi(x, y).dot(node2.delpsi(x, y));
        }
    }

    return diff1 * diff2 * ans;
}

long double MeshRect::integrateFPsi(Node node, long double a1, long double b1, long double a2, long double b2) {
    long double diff1 = (b1 - a1) / 2;
    long double avg1 = (b1 + a1) / 2;
    long double diff2 = (b2 - a2) / 2;
    long double avg2 = (b2 + a2) / 2;

    long double ans = 0;
    long double x,y;

    for (int i = 0; i < PRECISION; i++) {
        for (int j = 0; j < PRECISION; j++) {
            x = diff1 * abscissae[i] + avg1;
            y = diff2 * abscissae[j] + avg2;            
            ans += weights[i] * weights[j] * node.psi(x, y) * this->f(x, y);
        }
    }

    return diff1 * diff2 * ans;
}

void MeshRect::output(string fileName) {
    ofstream fout(fileName);

    fout << "# vtk DataFile Version 3.0" << endl;
    fout << "2D Poisson Equation Numeric Solution" << endl;
    fout << "ASCII" << endl;
    fout << endl;
    fout << "DATASET UNSTRUCTURED_GRID" << endl;
    fout << "POINTS " << Nb << " float" << endl;

    for (int n = 0; n < Nb; n++) {
        fout << (n % (Nx + 1)) * h1 << ' ' << (n / (Nx + 1)) * h2 << " 2" << endl;
    }

    fout << endl;
    fout << "CELLS " << Ne << " " << (Nlb + 1) * Ne << endl;

    for (int n = 0; n < Ne; n++) {
        fout << "4 ";

        for (int i = 0; i < Nlb; i++) {
            fout << Tb(i, n) << " ";
        }

        fout << endl;
    }

    fout << endl;
    fout << "CELL_TYPES " << Ne << endl;
    
    for (int n = 0; n < Ne; n++) {
        fout << 9 << endl;
    }

    fout << endl;
    fout << "POINT_DATA " << Nb << endl;
    fout << "SCALARS temprerature float 1" << endl;
    fout << "LOOKUP_TABLE default" << endl;

    for (int n = 0; n < Nb; n++) {
        fout << this->Uv(n) << endl;
    }

    fout.close();
}