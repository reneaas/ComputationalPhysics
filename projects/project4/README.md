# Project 4: The Ising model

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/report/Project_4.pdf).
The results were computed on three different computers:
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
To run the codes of the project, we advise to clone the repository pertaining to this project.

The structure of our codes are:
1. [main.cpp](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/main.cpp) contains the algorithms, generates data.
2. [main_mpi_MC.cpp](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/main_mpi_MC.cpp) contains the parallelized version of main.cpp.
3. [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/main.py) runs the cpp-files, write the produced data to files and categorize them into folders.
4. [plot.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/plot.py) read the files produced by main.py and generates plots from the data.


The codes are run in the following way:

## 1. Producing new results

The primary code of the numerical project is [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/main.py). It automates pretty much everything. A call to this can be done in the following way with respect to which problem one wants to solve:

The first thing you'll run if you wish to produce new data is

```console
python3 main.py
```

You'll then be prompted with:
```console
Which part of the project would you like to run? [b, c, e]
```

Where the letters correspond to:
- b : The 2x2 lattice case, comparison with analytical values.
- c : Estimation of equilibration & finding the probability distribution P(E)
- e : Phase transitions

Type in the letter corresponding to the part you wish to produce new data for. E.g if you want to produce new data for Estimation of equilibration, you'll type in

```console
c
```
### General for part b and c

Having chosen part "b" or "c" you'll be promted with the question

```console
Would you like an ordered or random initial spin matrix? [o/r]
```

Here you'll write "o" for an ordered initial spin matrix (ground state) or "r" for a randomized initial spin matrix.

Next you'll be promted with the question
```console
Specify number of Monte Carlo samples:
```

Here you type in the number of Monte Carlo samples you wish to run the code for. We suggest a number > 10^6, since the burn in period to reach equilibrium is by default t = 10^6. This period can be manually changed in main.cpp. In the report the burn in period is divided with total number of spins, corresponding to number of lattice sweeps.

### Part c

You'll now be prompted with the question
```console
Give temperature: [1.0 or 2.4]?
```
Here you write either "1.0" or "2.4" according to which temperature you wish to run the code for.

### Part e

Running part "e" will compile and run the parallelized Ising model. You can change the number of prosesses you would like to run in parallel by entering the [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/main.py) program and under part "e" changing the value of "p".

## 2. Graphic display of the results
All plotting of the produced data is done by running the program [plot.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/plot.py)

```console
python3 plot.py
```
You'll be promted with the question
```console
Which part of the project to run: [b, c, d, e, flags]
```

Where the letters correspond to:
##### b: The 2x2 lattice case

Plots expectation values for energy and magnetization, both with ordered and random initial spin matrices.

##### c : Estimation of equilibration

Choosing "c" you'll be prompted with the questions
```console
Give temperature: [1.0 or 2.4]
Number of Monte Carlo samples:
```

Type in the temperature you wish to plot results for, and the number of Monte Carlo samples used to produce the data you wish to plot.

##### d : Generating the probability distribution

As in part c you'll have to type in the temperature and number of Monte Carlo cycles.

Next you'll be prompted with the question
```console
ordered or random: [o/r]
```

Write either "o" or "r".

This now plots the probability distribution of energies, P(E), and the energy of the system as a function of Monte Carlo samples.


##### e : Phase transitions

This produces plots of the expectation values of energy, magnetization, heat capacity and magnetic susceptibility for multiple lattice sizes and temperatures. It also plots interpolation of the heat capacity for multiple lattice sizes and uses the maximum values to estimate the critical temperature, which is also plotted.



##### flags : Compiler flags

Plots a time comparison using different compiler flags.
