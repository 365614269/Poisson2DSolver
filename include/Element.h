#pragma once
#include "Node.h"
#include "globals.h"
#include <vector>

using std::vector;
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