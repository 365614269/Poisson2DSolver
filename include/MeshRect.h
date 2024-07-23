#pragma once

#include "Element.h"
#include "abscissae.h"
#include "weights.h"
#include "exprtk.hpp"
#include "rapidxml.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include "globals.h"
#include "rapidxml_utils.hpp"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef exprtk::symbol_table<double> symbol_table;
typedef exprtk::expression<double>   expression;
typedef exprtk::parser<double>       parser;

class MeshRect {
    private : 
        vector<Element> elements;
        MatrixXd stiffness;
        VectorXd Uv,delUv,Fv;
        expression f_expr, delf_expr;

        double u_val;

        expression parse_f();
        expression parse_delf();

        double u(double, double);
        double f(double, double);
        double delf(double, double);
        double g(double, double);

        static vector<string> split(string, string);

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