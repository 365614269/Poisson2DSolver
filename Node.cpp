#include "Node.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

Node::Node(long double h_1, long double h_2, long double N_y, int i_, int j_, long double (*f_)(long double, long double), long double (*delf_)(long double, long double)) {
    this->h1 = h_1;
    this->Ny = N_y;
    this->h2 = h_2;
    this->i = i_;
    this->j = j_;
    this->f = f_;
    this->delf = delf_;
}

long double Node::psi(long double x, long double y, int localNodeIndex) {
    int sign[4][2] = {{-1, 1}, {1, 1}, {1, -1}, {-1, -1}};

    int signs = sign[localNodeIndex][0];
    int signt = sign[localNodeIndex][1];

    long double x0 = (this->j + 0.5) * this->h1;
    long double y0 = (this->Ny * this->h2) - (this->i + 0.5) * this->h2;
    long double s = (2.0 / this->h1) * (x - x0);
    long double t = (2.0 / this->h2) * (y - y0);

    return 0.25 * (1 + signs * s) * (1 + signt * t);
}

VectorXd Node::delpsi(long double x, long double y, int localNodeIndex) {
    int sign[4][4] = {{-1, -1, 1, -1}, 
                      {1, 1, 1, 1},
                      {1, -1, -1, -1},
                      {-1, 1, -1, 1}};

    VectorXd v(2);

    long double x0 = (this->j + 0.5) * this->h1;
    long double y0 = (this->Ny * this->h2)- (this->i + 0.5) * this->h2;
    long double s = (2.0 / this->h1) * (x - x0);
    long double t = (2.0 / this->h2) * (y - y0);

    int sign1 = sign[localNodeIndex][0];
    int sign2 = sign[localNodeIndex][1];
    int sign3 = sign[localNodeIndex][2];
    int sign4 = sign[localNodeIndex][3];

    v << ((sign1 + sign2 * t) / (2 * this->h1)), ((sign3 + sign4 * s) / (2 * this->h2));

    return v;
}

long double Node::int1(long double x, long double y, int alpha, int beta) {
    return this->delpsi(x, y, alpha).dot(this->delpsi(x, y, beta));
}

long double Node::int2(long double x, long double y, int alpha, int beta) {
    return this->psi(x, y, alpha) * this->psi(x, y, beta) * delf(x, y);
}

long double Node::int3(long double x, long double y, int alpha, int beta) {
    return this->psi(x, y, alpha) * this->f(x, y);
}