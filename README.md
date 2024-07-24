# Introduction

The project aims to use Finite Element Method (FEM) to find numeric solutions for the 2D non-linear Poisson equation defined in the rectangle $$(0, 0)$$ to $$(lx, ly)$$:

$$\nabla\cdot\nabla u = f(u)$$

i.e.

$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = f(u)$$

where $$u$$ is a 2D function for us to find, and $$f$$ is a function of $$u$$.

Note that $$u = g(x, y)$$ on the boundary as Dirichlet Boundary Condition and $$u = u(x, y)$$  elsewhere.

Neumann Boundary Condition is not considered in this project.

# Getting started

No dependencies other than your C++ compiler is needed. 

To compile the code, go to the cloned directory and compile every .cpp file in the directory.

The command may vary for different compilers. Here's the command for g++:

`$ g++ -o main *.cpp`

However, if you use a windows or mac operating system, there are releases.

When you run the executable, give the directory of the configuration file as a command-line argument.

`$ ./your_executable /path/to/your_config.xml`

If you just want to use the template file **data.xml**, use the command 

`$ ./your_executable datas.xml`

# Configuration

At this point, you might be curious about how to configure the program to let it solve the equation you want. However, if you just want to play around with the program, skip this part.

$$lx$$ and $$ly$$ represents the rectangular area where the equation is defined.

$$Nx$$ and $$Ny$$ represents the number of elements per row and column.

The boundary nodes and their values are given in two strings respectively, acting as Dirichlet Boundary Condition. If a boundary node's value is not specified, it will be considered as a free node.

The initial guess function $$u_0(x, y)$$, source function $$f(u)$$ and derivative of source function $$\frac{\partial f}{\partial u}$$ are represented as a string and parsed to the program.

The absolute tolerance of error and relative tolerance of error are specified.

The directory of output file is given at the end.

# Visualize the result

The result must be stored as a vtk file, so that the output file directory must end with ".vtk".

Download [Paraview](https://www.paraview.org/download/) to draw a diagram of the result. Simply load the output vtk file, and apply it. You will see a colored rectangle representing the result.

# End
