# Project 4: The Ising model

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project3/report/Project3_report.pdf).
The results were computed on three different computes.:
1. Computer 1:
  - CPU: Intel i7-8565U (8)
  - RAM: 16 GB

2. Computer 2:
  - CPU: Intel i5-7360U (4)
  - RAM: 8 GB

3. Computer 3:
  - CPU: Intel i5-4260U (4)
  - RAM: 4 GB


## Codes (Documentation)
To run the codes of the project, we advise to clone the repository pertaining to this project and run the codes in the following way:

The primary code of the numerical project is [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project3/codes/main.py). It automates pretty much everything. A call to this can be done in the following way with respect to which problem one wants to solve:

The first thing you'll run if you wish to produce new data is

```console
python3 main.py
```

* Once this particular problem is run, you'll be prompted with:
```console
Which part of the project would you like to run? [b, c, d, e]
```

Where the letters correspond to:
- b : 2x2-lattice
- c : Estimation of equilibration
- d : The probability distribution P(E)
- e : Phase transitions

Type in the letter corresponding to the part you wish to produce new data for. E.g if you want to produce new data for Estimation of equilibration, you'll type in

```console
c
```


b. 2x2-lattice

Having chosen part "b", you'll be prompted with the question

```console
Would you like an ordered or random initial spin matrix? [o/r]
```

Here you'll write "o" for an ordered initial spin matrix (ground state) or "r" for a randomized initial spin matrix.

Next you'll be promted with the question
```console
Specify number of Monte Carlo samples:
```

Here you



3. Benchmarking the Laguerre method
```console
python3 main.py benchmark_laguerre
```
* Once this particular problem is run, you'll be prompted with:
```console
Produce new data? Type yes or no:
```
  - Type yes to compile and produce new results and new plots. Type no to not compile and to produce results and plots with data already saved.

4. Comparing the Legendre and the Laguerre methods
```console
python3 main.py compare_gauss
```

* Running this produces plots that compares the relative error and CPU times of the gaussian quadrature methods.

5. Benchmarking the Monte Carlo methods and Ground State energy
```console
python3 main.py multiple_MC
```

* By running this command you will compile the Brute force, the Brute force with MPI, the Importance sampling and the Importance sampling with MPI versions of the Monte Carlo method. You'll then be prompted with:
```console
Produce new data? Type yes or no:
```

  - Type yes to compile and produce new sets of data. Type no to use already existing datasets for sample sets N = [100,...,10000] with a stepsize dN = 100.
  - If you type yes, you'll be prompted with
    ```console
    How many sets of data do you want?:
    ```
    - Here you choose how many times you want to run and produce datasets
  - If you type no, you print results for m = 10 000 datasets (which is the current data stored in the repository).

6. Estimating infinity (Lambda)
```console
python3 main.py find_lambda
```

* Once this particular problem is run, you'll be prompted with:
```console
Produce new data? Type yes or no:
```
 * If you choose yes you'll compile the program and run the Laguerre method for several lambdas to estimate which value of lambda gives the lowest relative error for integration points N = 25 (If you wish, you can change N in the main.py program under "find_lambda"). If you choose no, you'll get the lambda which gives the lowest relative error for results already produced for N = 25.
