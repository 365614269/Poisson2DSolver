# Introduction

The project aims to use Finite Element Method (FEM) to find numeric solutions for the 2D non-linear Poisson equation defined in the rectangle $$(0, 0)$$ to $$(lx, ly)$$ defined in the 2D $$(x, y)$$ plane:

$$\nabla\cdot\nabla u = f(u)$$

i.e.

$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = f(u)$$

where $$u$$ is a 2D function for us to find, and $$f$$ is a function of $$u$$.

Note that $$u = g(x, y)$$ on the boundary as Dirichlet Boundary Condition and $$u = u(x, y)$$  elsewhere.

Neumann Boundary Condition is not considered in this project.

# Getting started

No dependencies other than your C++ compiler is needed. 

To compile the code, go to the cloned directory and compile every .cpp file in the directory.

However, if you use a windows or mac operating system, there are releases.

The command may vary for different compilers. Here's the command for g++:

`$ g++ -O3 -Wa,-mbig-obj -o your_executable *.cpp`

Since the program requires intense calculation, you should always compile with the highest optimization level.

When you run the executable, give the directory of the configuration file as a command-line argument.

`$ ./your_executable /path/to/your_config.xml`

If you just want to use the template file **datas.xml**, use the command 

`$ ./your_executable datas.xml`

# Configuration

At this point, you might be curious about how to configure the program to let it solve the equation you want. However, if you just want to play around with the program, skip this part.

$$lx$$ and $$ly$$ represents the rectangular area where the equation is defined.

$$Nx$$ and $$Ny$$ represents the number of elements per row and column.

The boundary conditions are specified through four lines, with line name="u_XY" and value="AAA-BBB,VAL1;CCC-DDD,VAL2". The first part "AAA-BBB,VAL1" representing the interval of the percentage AAA-BBB of nodes from X to Y set to VAL1.

The initial guess function $$u_0(x, y)$$, source function $$f(u)$$ and derivative of source function $$\frac{\partial f}{\partial u}$$ are represented as a string and parsed to the program.

The absolute tolerance of error and relative tolerance of error are specified.

The directory of output file is given at the end.

# Visualize the result

The result must be stored as a vtk file, so that the output file directory must end with ".vtk".

Download [ParaView](https://www.paraview.org/download/) to draw a diagram of the result. Simply load the output vtk file, and apply it. You will see a colored rectangle representing the result.

# End
