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
Here you compile the code by typing "no_mpi". If you have already compiled it, and just want to change some of the input parameters, you call:
```console
python3 main.py 0
```

  * Once this particular problem is run (either with no_mpi or 0), you'll be prompted with:
  ```console
Specify integration method, choose from:
Gauss Legendre method                  --> type 1
Gauss Laguerre method                  --> type 2
Brute force Monte Carlo                --> type 3
Monte Carlo with importance sampling   --> type 4
```
- If you type 1, you'll have to specify number of integration points (an integer), and then the integration limits [a,b] (two double floating point numbers).
- If you type 2, you will have to specify number of integration points (an integer). Then you'll have to decide if you want to do the integral for three or six dimensions:
```console
Choose dimension for integral:
For 3 dimensions                         --> type 3
For 6 dimensions                         --> type 6
```
You then obtain the exact and calculated value for the integral.
- If you type 3, you'll have to specify number of monte carlo samples (an integer), then the integration limits [a,b] (two double floating point numbers). The results appear in the terminal.
- If you type 4, you'll have to specify number of monte carlo samples (an integer), then the maximum radial distance (double floating point number). The results appear in the terminal.

2. Benchmarking the Laguerre method
```console
python3 main.py benchmark_laguerre
```
* Once this particular problem is run, you'll be prompted with:
```console
Produce new data? Type yes or no:
```
  - Type yes to compile and produce new results and new plots. Type no to not compile and to produce results and plots with data already saved.

3. Benchmarking the Monte Carlo methods and Ground State energy
```console
python3 main.py multiple_MC
```

* By running this command you will compile the Brute force, the Brute force with MPI, the Importance sampling and the Importance sampling with MPI versions of the Monte Carlo method. You'll then be prompted with:
```console
Produce new data? Type yes or no:
```

  - Type yes to compile and produce new sets of data. Type no to use already existing datasets for sample sets N = [10,10^2,10^3,10^4,10^5,10^6] with m = 10 000 datasets (Don't try this at home).
  - If you type yes, you'll be prompted with
    ```console
    How many sets of data do you want?:
    ```
    - Here you choose how many times you want to run and produce datasets for   sample sets  N = [10,10^2,10^3,10^4,10^5,10^6].
  - If you type no, you print results for m = 10 0000 datasets.
