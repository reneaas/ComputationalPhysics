# Project 1: Tridiagonal matrix solver

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project5/report/Project_5.pdf)
outlines the mathemathical background of the diffusion equation.


## Codes (Documentation)
The [C++ code](./codes/cpp) included in this project is solving the 2D diffusion equation using an explicit scheme.

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
