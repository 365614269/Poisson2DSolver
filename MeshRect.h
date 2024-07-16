#pragma once
#include <iostream>
#include <utility>
#include "eigen/Eigen/Dense"
#include "integrate2D.h"
#include "Node.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class MeshRect {
    private : 
        long double lx,ly;
        long double (*f) (long double, long double);
        long double (*delf) (long double, long double);
        int Nx,Ny,Nb;
        Node*** elements;
        MatrixXd stiffness;
        VectorXd Uv,Fv;
    public :
        MeshRect(long double, long double, int, int, long double (*)(long double, long double), long double (*)(long double, long double), VectorXd);
        // long double u(long double, long double);
        // pair<int, int> elem(long double, long double);
        Node getNode(int, int);
        int exchangeIndex(int, int);
        pair<int, int> exchangeIndex(int);
        MatrixXd getStiffness();
        VectorXd getF();
        VectorXd getU();
        void setU(VectorXd);
        int Tb(int, int, int);
        void addAF(Node);
        void calculateAF();
        void displayStiffness();
        void displayF();
        void applyBC();
};