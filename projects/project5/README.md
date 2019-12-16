# Project 5: The Diffusion Equation

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project5/report/Project_5.pdf).
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
2. [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/main.py) runs the cpp-files, write the produced data to files and categorize them into folders.
3. [plot.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/plot.py) read the files produced by main.py and generates plots from the data.


The codes are run in the following way:

## 1. Producing new results

The primary code of the numerical project is [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project4/codes/main.py). It automates pretty much everything. A call to this can be done in the following way with respect to which problem one wants to solve:

The first thing you'll run if you wish to produce new data is

```console
python3 main.py
```

You'll then be prompted with:
```console
__________________________________________
Which part would you like to run?
1-dimensional schemes          --> Type 1
2-dimensional scheme           --> Type 2
Stability analysis for 1D-case --> Type 3
__________________________________________
```

Type in the number corresponding to the part you wish to produce new data for. E.g if you want to produce new data for the 1-dimensional schemes you'll type in

```console
1
```

### 1-dimensional schemes
Choosing "1" produces results for dx = 0.1 and dx = 0.01 and simulation time t = 0.02 and t = 1, for r = Δt/Δx^2 = 0.5.

### 2-dimensional scheme
Choosing "2" produces results for steplength h = 0.01, simulation time t = 0.5 and r = Δt/Δx^2 = 0.25.

### Stability analysis for 1D-case
Choosing "3" produces a stability analysis of the three 1D-schemes comparing results for r = 0.5 and r = 0.505.

## 2. Plotting results
The first thing you'll run if you wish to plot data is

```console
python3 plot.py
```

You'll then be prompted with:
```console
______________________________________________________
Which part would you like to plot?
1-dimensional schemes                    --> Type 1
2-dimensional scheme                     --> Type 2
Stability analysis for 1D-schemes        --> Type 3
Contour plot for Crank-Nicolson scheme   --> Type 4
______________________________________________________

```

Type the number corresponding to the results you wish to produce plots for. E.g if you want to produce plots for the 1-dimensional schemes you'll type in

```console
1
```
