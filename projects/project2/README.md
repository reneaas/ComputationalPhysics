# Project 2: Eigenvalue problems

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/report/Project_2_report.pdf)
provides the necessary mathematical background to understand the Jacobi rotation algorithm and showcases application
of this algorithm to two second-order differential equations.

## Code documentation
The [C++ code](./codes/cpp) included in this project solves a specified eigenvalue equation
using the Jacobi rotation algorithm.

### Dependencies

The code project itself utilizes [OpenMP](https://www.openmp.org) and [Armadillo](http://arma.sourceforge.net).
To build the C++ codes, then, these libraries must be installed.

### Building the code

The [C++ codes](./codes/cpp) can be built by

```sh
make compile && make link
```

To execute the code, you can utilize the makefile by

```sh
make run
```

To build and execute all in one go, execute the following in the command line:

```sh
make all
```
