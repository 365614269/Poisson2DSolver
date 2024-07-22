#pragma once

#include "eigen/Eigen/Dense"
#include <utility>

using Eigen::VectorXd;
using namespace std;

#define PRECISION 10
#define MAX_ITERATIONS 10
#define ABSOLUTE_TOLERANCE 10e-4
#define RELATIVE_TOLERANCE 10e-4

const double lx = 2;
const double ly = 1;
const int Nx = 20;
const int Ny = 10;
const int Nb = (Nx + 1) * (Ny + 1);
const int Ne = Nx * Ny;
const int Nlb = 4;
const double h1 = lx / Nx;
const double h2 = ly / Ny;
const string shape = "rectangle";