#include "MeshRect.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

MeshRect::MeshRect(long double l_x, long double l_y, int N_x, int N_y, long double (*f_)(long double, long double), long double (*delf_)(long double, long double), VectorXd U_) {
    this->lx = l_x;
    this->ly = l_y;
    this->Nx = N_x;
    this->Ny = N_y;
    this->f = f_;
    this->delf = delf_;
    this->Uv = U_;
    this->Fv = VectorXd::Zero((Nx + 1) * (Ny + 1));
    this->elements = new Node**[Ny];
    this->stiffness = MatrixXd::Zero((this->Nx + 1) * (this->Ny + 1), (this->Nx + 1) * (this->Ny + 1));

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

void MeshRect::addStiffness(Node node) {
    int Nlb = 4;
    long double ans1, ans2, ans3;

    for (int alpha = 0; alpha < Nlb; alpha++) {
        for (int beta = 0; beta < Nlb; beta++) {
            long double x1,x2,y1,y2;

            x1 = node.j * node.h1;
            x2 = (node.j + 1) * node.h1;
            y1 = (node.Ny - node.i - 1) * node.h2;
            y2 = (node.Ny - node.i) * node.h2;
            ans1 = integral(&Node::int1, alpha, beta, x1, x2, y1, y2, node);
            ans2 = integral(&Node::int2, alpha, beta, x1, x2, y1, y2, node);

            long double ans = -ans1-ans2;

            this->stiffness(Tb(beta, node.i, node.j), Tb(alpha, node.i, node.j)) += ans;
        }
    }
}

void MeshRect::calculateStiffness() {
    for (int i = 0; i < this->Ny; i++) {
        for (int j = 0; j < this->Nx; j++) {
            Node node = this->getNode(i, j);
            this->addStiffness(node);
        }
    }
}


void MeshRect::displayStiffness() {
    cout << this->stiffness << endl;
}

int MeshRect::exchangeIndex(int x, int y) {
    return x * this->Nx + y;
}