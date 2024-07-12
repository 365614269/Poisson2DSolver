#pragma once

#include "eigen/Eigen/Dense"

using namespace std;


class Node {
    public :
        long double h1, h2;
        int Ny;
        int i,j;
        long double (*delf)(long double, long double);

        Node (long double, long double, long double, int i, int j, long double (*)(long double, long double));
        long double psi(long double, long double, int);
        Eigen::VectorXd delpsi(long double, long double, int);
        long double int1(long double, long double, int, int);
        long double int2(long double, long double, int, int);
};