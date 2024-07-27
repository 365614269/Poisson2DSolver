#pragma once

#include "globals.h"

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