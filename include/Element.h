#pragma once
#include "Node.h"

using Eigen::VectorXd;

class Element {
    private :
        int elemIndex;
        Node nodes[Nlb];

    public :
        Element(int);
        Element(); // Default constructor
        int getElemIndex();
        Node& getNode(int);
};