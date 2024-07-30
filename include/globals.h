#pragma once

#include "../eigen/Dense"
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

extern string u_AB;
extern string u_BC;
extern string u_DC;
extern string u_AD;