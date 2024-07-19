#include "Node.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Node::Node(int elemIndex, int localNodeIndex) {
    this->elemIndex = elemIndex;
    this->localNodeIndex = localNodeIndex;
}

Node::Node() {}


long double Node::psi(long double x, long double y) {
    int sign[4][2] = {{-1, 1}, {1, 1}, {1, -1}, {-1, -1}};

    int signs = sign[this->localNodeIndex][0];
    int signt = sign[this->localNodeIndex][1];

    pair<int, int> gridIndex = make_pair(this->elemIndex / Nx, this->elemIndex % Nx);  // exchange index

    long double x0 = (gridIndex.second + 0.5) * h1;
    long double y0 = (Ny * h2) - (gridIndex.first + 0.5) * h2;
    long double s = 2 * (x - x0) / h1;
    long double t = 2 * (y - y0) / h2;

    return 0.25 * (1 + signs * s) * (1 + signt * t);
}

VectorXd Node::delpsi(long double x, long double y) {
    int sign[4][4] = {{-1, -1, 1, -1}, 
                      {1, 1, 1, 1},
                      {1, -1, -1, -1},
                      {-1, 1, -1, 1}};

    VectorXd v(2);

    pair<int, int> gridIndex = make_pair(this->elemIndex / Nx, this->elemIndex % Nx);  // exchange index

    long double x0 = (gridIndex.second + 0.5) * h1;
    long double y0 = (Ny * h2) - (gridIndex.first + 0.5) * h2;
    long double s = 2 * (x - x0) / h1;
    long double t = 2 * (y - y0) / h2;

    int sign1 = sign[this->localNodeIndex][0];
    int sign2 = sign[this->localNodeIndex][1];
    int sign3 = sign[this->localNodeIndex][2];
    int sign4 = sign[this->localNodeIndex][3];

    v << ((sign1 + sign2 * t) / (2 * h1)), ((sign3 + sign4 * s) / (2 * h2));

    return v;
}