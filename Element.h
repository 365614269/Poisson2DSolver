#pragma once
#include "Node.h"

using Eigen::VectorXd;

class Element {
    public :
        int elemIndex;
        Node nodes[Nlb];

        Element(int);
        Element(); // Default constructor
};