# Project 3: Numerical integration

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/report/Project_2_report.pdf). SETT INN LINK TIL RIKTIG RAPPORT

## Codes (Documentation)
To run the codes of the project, I advise to clone the repository pertaining to this project and run the codes in the following way:

The primary code of the numerical project is [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project3/codes/main.py). It automates pretty much everything. A call to this can be done in the following way with respect to which problem one wants to solve:

1. Gauss-Legendre-, Gauss-Laguerre quadrature, brute force Monte Carlo or improved Monte Carlo without parallization:

```console
python3 main.py no_mpi
```

  * Once this particular problem is run, you'll be prompted with:
    1. Specify integration method, choose from:
        ------------------------------------------------------
        Gauss Legendre method                  --> type 1
        Gauss Laguerre method                  --> type 2
        Brute force Monte Carlo                --> type 3
        Monte Carlo with importance sampling   --> type 4
        ------------------------------------------------------

        - If you type 1, you'll have to specify number of integration points (an integer), and then the integration limits [a,b] (two double floating point numbers).
        - If you type 2, you will have to specify number of integration points (an integer). Then you'll have to decide if you want to do the integral for three or six dimensions:

        Choose dimension for integral:
        ------------------------------------------------------
        For 3 dimensions                         --> type 3
        For 6 dimensions                         --> type 6
        ------------------------------------------------------

        You then obtain the exact and calculated value for the integral.
        - If you type 3, you'll have to specify number of monte carlo samples (an integer), then the integration limits [a,b] (two double floating point numbers). The results appear in the terminal.
        - If you type 4, you'll have to specify number of monte carlo samples (an integer), then the maximum radial distance (double floating point number). The results appear in the terminal.
