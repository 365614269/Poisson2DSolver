#pragma once

#include "../Eigen/Dense"
#include <utility>

using Eigen::VectorXd;
using namespace std;

#define PRECISION 10
#define MAX_ITERATIONS 10

extern double lx;
extern double ly;
extern int Nx;
extern int Ny;
extern string shape;

extern int Nb;
extern int Ne;
extern int Nlb;
extern double h1;
extern double h2; 

extern double abs_tol;
extern double rel_tol;

extern string u0_str;
extern string f_str;
extern string delf_str;
extern string output_file_dir;

extern string BCValues;
extern string BCNodes;