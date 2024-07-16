#include "MeshRect.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

MeshRect::MeshRect(long double l_x, long double l_y, int N_x, int N_y, long double (*f_)(long double, long double), long double (*delf_)(long double, long double), VectorXd U_) {
    this->lx = l_x;
    this->ly = l_y;
    this->Nx = N_x;
    this->Ny = N_y;
    this->Nb = (N_x + 1) * (N_y + 1);
    this->f = f_;
    this->delf = delf_;
    this->Uv = U_;
    this->Fv = VectorXd::Zero(Nb);
    this->elements = new Node**[Ny];
    this->stiffness = MatrixXd::Zero(Nb, Nb);

    for(int i = 0; i < Ny; ++i){
        this->elements[i] = new Node*[this->Nx];

        for(int j = 0; j < Nx; ++j){
            this->elements[i][j] = new Node(this->lx / this->Nx, this->ly / this->Ny, this->Ny, i, j, this->f, this->delf);
        }
    }
}

Node MeshRect::getNode(int i, int j) {
    return *(this->elements[i][j]);
}

int MeshRect::Tb(int localNodeIndex, int i, int j) {
    int offset[4][2] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
    int x = i + offset[localNodeIndex][0];
    int y = j + offset[localNodeIndex][1];

    return x * (this->Nx + 1) + y;
}

void MeshRect::addAF(Node node) {
    int Nlb = 4;
    long double ans1, ans2, ans3;
    long double x1,x2,y1,y2;

    x1 = node.j * node.h1;
    x2 = (node.j + 1) * node.h1;
    y1 = (node.Ny - node.i - 1) * node.h2;
    y2 = (node.Ny - node.i) * node.h2;

    for (int alpha = 0; alpha < Nlb; alpha++) {
        ans3 = 0;

        for (int beta = 0; beta < Nlb; beta++) {
            ans1 = integral(&Node::int1, alpha, beta, x1, x2, y1, y2, node);
            ans2 = integral(&Node::int2, alpha, beta, x1, x2, y1, y2, node);
            ans3 += ans1 * this->Uv(Tb(alpha, node.i, node.j));

            long double ans = -ans1-ans2;

            this->stiffness(Tb(beta, node.i, node.j), Tb(alpha, node.i, node.j)) += ans;
        }
        
        ans3 += integral(&Node::int3, alpha, alpha, x1, x2, y1, y2, node);
        this->Fv(Tb(alpha, node.i, node.j)) += ans3;
    }
}

void MeshRect::calculateAF() {
    for (int i = 0; i < this->Ny; i++) {
        for (int j = 0; j < this->Nx; j++) {
            Node node = this->getNode(i, j);
            this->addAF(node);
        }
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

int MeshRect::exchangeIndex(int x, int y) {
    return x * this->Nx + y;
}

pair<int, int> MeshRect::exchangeIndex(int n) {
    return make_pair(n / (this->Nx + 1), n % (this->Nx + 1));
}

void MeshRect::applyBC() {
    for (int n = 0; n < this->Nb; n++) {
        pair<int, int> index = this->exchangeIndex(n);

        if (index.first == 0 || index.first == this->Ny || index.second == 0 && index.second == this->Nx) {
            for (int i = 0; i < this->Nb; i++) {
                this->stiffness(n,i) = 0;
            }
            this->stiffness(n,n) = 1;
            this->Fv(n) = 0;
        }
    }
}

// pair<int, int> MeshRect::elem(long double x, long double y) {
//     long double h1 = this->lx / this->Nx;
//     long double h2 = this->ly / this->Ny;

//     int a = (this->ly - y) / h2;
//     int b = x / h1;

//     a = Math.max(this->Ny, a);
//     b = Math.max(this->Nx, b);

//     return make_pair(a, b);
// }

// long double MeshRect::u(long double x, long double y){
//     //find which element (x, y)is in.
//     pair<int, int> index = this->elem(x, y);
//     Node node = this->getNode(index.first, index.second);

//     int Nb = 4;
//     VectorXd v1(Nb);
//     VectorXd v2(Nb);

//     // add up the 4 components * shape functions of the 4 nodes
//     for (int i = 0; i < Nb; i++) {
//         v1 << node.psi(x, y, Nb);

//         int n = Tb(i, index.first, index.second);
//         v2 << this->Uv(n);
//     }

//     return v1.dot(v2);
// }