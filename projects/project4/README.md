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
Which part of the project would you like to run? [b, c, e]
```

Where the letters correspond to:
- b : 2x2-lattice
- c : Estimation of equilibration & The probability distribution P(E)
- e : Phase transitions

Type in the letter corresponding to the part you wish to produce new data for. E.g if you want to produce new data for Estimation of equilibration, you'll type in

```console
c
```
### General for b and c

Having chosen part "b", "c" or "d" you'll be promted with the question

```console
Would you like an ordered or random initial spin matrix? [o/r]
```

Here you'll write "o" for an ordered initial spin matrix (ground state) or "r" for a randomized initial spin matrix.

Next you'll be promted with the question
```console
Specify number of Monte Carlo samples:
```

Here you type in the number of Monte Carlo samples you wish to run the code for. We suggest number > 10^6.

### Part c

You'll now be promted with the question
```console
Run for temperature 1 or 2.4?
```
Here you write either "1" or "2.4" according to which temperature you wish to run the code for.

### Part e

Running part "e" will compile and run the parallelized Ising model. 
