# Project 2: Eigenvalue problems

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/report/project2.pdf).

## Codes (Documentation)
To run the codes of the project, I advise to clone the repository pertaining to this project and run the codes in the following way:

The primary code of the numerical project is [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/codes/main.py). It automates pretty much everything. A call to this can be done in the following way with respect to which problem one wants to solve:

1. Buckling Beam problem:

```console
python3 main.py Number_of_gridpoints Max_number_of_iterations bb compile_instruction
```

2. Quantum dots with one electron:

```console
python3 main.py Number_of_gridpoints Max_number_of_iterations qm1 compile_instruction
```

  * Once this particular problem is run, you'll be promted with:
    1. Do you wish to run this for several infinities? Type yes or no:
      - If you type yes, you'll run the code for several predefined "infinities", which in this case is to be interpreted as the maximum value of rho.
      - if you type no, you'll be prompted to manually insert a particular value of rho which you wish to run the simulation with.

3. Quantum dots with two electrons:

```console
python3 main.py Number_of_gridpoints Max_number_of_iterations qm2 compile_instruction
```

  * Here you will be prompted to insert additional values for :
    1. Repulsion: Whether to include electron repulsion or not.
    2. Angular frequency: A double floating point number.


By compilation_instruction, you simply type "yes" if you wish to recompile the code. If recompilation is undesired, type "no".


As an example, say you want to solve the quantum mechanics problem involving two electrons with repulsion with the following parameters:
* Number_of_gridpoints = 100
* Max_number_of_iterations = 100000

Furthermore assume that you want to compile the code, then the code can be run in the following way:

```console
python3 main.py 100 100000 qm2 yes
```
The program will now prompt you with several inputs. The example below shows how we to correctly respond if repulsion is to be included and the angular
frequency is to be set to 0.5:
```console
Solving Schrödingers eq in 3D with two electrons.
Running code for n = 100
Include electron repulsion? Type yes or no: yes
Give the angular frequency: 0.5
```

Here we give a summary of the codes:
- [main.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/codes/main.py): main code that runs the necessary codes to perform computations, and organizes and moves files around.
- [main.cpp](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/codes/main.cpp): is the main program where computations happen. It facilitates the use of functions from functions.cpp to run the main algorithms.
- [functions.cpp](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/codes/functions.cpp): contains all functions written to make the main program run.
- [functions.h](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/codes/functions.h): Header file
- [make_plot.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/codes/make_plot.py): Creates a plot of the ground state wavefunction based on the computations done in main.cpp
- [plot_wavefunctions.py](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/codes/plot_wavefunctions.py): plots several wavefunctions together where the only parameter that differs is the angular frequency
