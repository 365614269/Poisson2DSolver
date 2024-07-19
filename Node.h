#pragma once

#include <iostream>
#include "eigen/Eigen/Dense"
#include "globals.h"

using namespace std;
using Eigen::VectorXd;

class Node {
    public :
        int elemIndex;
        int localNodeIndex;
        Node(int, int);
        Node(); // Default constructor
        long double psi(long double x, long double y);
        VectorXd delpsi(long double x, long double y);
};