#pragma once
#include "Node.h"
#include "globals.h"

using Eigen::VectorXd;

class Element {
    private :
        int elemIndex;
        vector<Node> nodes;

    public :
        Element(int);
        Element(); // Default constructor
        int getElemIndex();
        Node& getNode(int);
};