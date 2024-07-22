#pragma once

#include <iostream>
#include "eigen/Eigen/Dense"
#include "globals.h"

using namespace std;
using Eigen::VectorXd;

class Node {
    private : 
        int elemIndex;
        int localNodeIndex;

    public :
        Node(int, int);
        Node(); // Default constructor
        double psi(double x, double y);
        VectorXd delpsi(double x, double y);
};