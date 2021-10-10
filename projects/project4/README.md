# Project 4: The Ising model

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/report/Project_4.pdf)
outlines the necessary mathematical background to implement Monte Carlo simulations of the 2D Ising model
using the Metropolis-Hastings algorithm to generate samples. The endgame of the project was to
determine the critical temperature of the 2D Ising model numerically, which you can read more about in the report.


## Codes

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
