# Project 1: Tridiagonal matrix equations

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project5/report/Project_1.pdf)
provides a review of the mathematics behind the Thomas Algorithm,
an algorithm designed to solve matrix equations Ax = b where A is tridiagonal.
We also provide a review of LU-decomposition specialized for tridiagonal matrices.


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
