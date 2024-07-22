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

        double u(double, double);
        double f(double, double);
        double delf(double, double);
        double g(double, double);

        pair<int, int> elem(double, double);

        Element& getElement(int);

        int Tb(int, int, int);
        int Tb(int, int);

        double integratePsiDelF(Node&, Node&, double, double, double, double);
        double integrateDelpsi(Node&, Node&, double, double, double, double);
        double integrateFPsi(Node&, double, double, double, double);

    public :
        MeshRect(VectorXd&);
        void applyBCtoU();
        void applyBCtoDelU();
        void addAF(int);
        void calculateAF();
        void addDelU();

        bool endConditionMet();
        void output(string);
};