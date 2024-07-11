#pragma once
#include <iostream>
#include <utility>
#include "eigen/Eigen/Dense"
#include "integrate2D.h"
#include "Node.h"

class MeshRect {
    private : 
        long double lx,ly;
        // long double (*f) (long double, long double);
        // long double (*delf) (long double, long double);
        int Nx,Ny;
        Node*** elements;
        Eigen::MatrixXd stiffness;
    public :
        MeshRect(long double, long double, int, int);
        Node getNode(int, int);
        int exchangeIndex(int, int);
        int Tb(int, int, int);
        void addStiffness(Node);
        void calculateStiffness();
        void displayStiffness();
};