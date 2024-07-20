#pragma once

#include "Element.h"
#include "abscissae.h"
#include "weights.h"
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class MeshRect {
    private : 
        Element elements[Ne];
        MatrixXd stiffness;
        VectorXd Uv,delUv,Fv;

    public :
        MeshRect(VectorXd);
        long double u(long double, long double);
        long double f(long double, long double);
        long double delf(long double, long double);
        long double g(long double, long double);

        pair<int, int> elem(long double, long double);

        Element getElement(int);

        int Tb(int, int, int);
        int Tb(int, int);

        long double integratePsiDelF(Node, Node, long double, long double, long double, long double);
        long double integrateDelpsi(Node, Node, long double, long double, long double, long double);
        long double integrateFPsi(Node, long double, long double, long double, long double);

        void applyBCtoU();
        void applyBCtoDelU();
        void addAF(int);
        void calculateAF();
        void addDelU();

        bool endConditionMet();
        void output(string);
};