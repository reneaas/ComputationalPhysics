# Project 2: Eigenvalue problems

## The report
The [report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/report/project2.pdf).

## Codes (Documentation)
To run the codes of the project, I advise to clone the repository pertaining to this project and run the codes in the following way:

### main.py
This is the main code of the project. It essentially automates everything. What you need to write in a standard Unix console is:

1. Buckling Beam problem:

```console
python3 main.py Number_of_gridpoints Max_number_of_iterations bb
```

2. Quantum dots with one electron:

```console
python3 main.py Number_of_gridpoints Max_number_of_iterations qm1
```

3. Quantum dots with two electrons:

```console
python3 main.py Number_of_gridpoints Max_number_of_iterations qm2
```

  * Here you will be prompted to insert additional values for :
    1. Repulsion: Whether to include electron repulsion or not.
    2. Angular frequency: A double floating point number.

On every run, you'll be prompted with whether you want to compile the code again. If no changes has been made to the code, type "no". Otherwise type "yes".

As an example, say you want to solve the quantum mechanics problem involving two electrons with repulsion and let's say we use the following parameters:
* Number_of_gridpoints = 100
* Max_number_of_iterations = 100000
Then you can run the code as

```console
python3 main.py 100 100000 qm2
```
Now assume we do not want to compile the code, we want to include repulsion and run the code with angular frequency set to 0.5. Then the console would look like this:

```console
Compile anew? Type yes or no: no
Solving Schrödingers eq in 3D with two electrons.
Running code for n = 100
Include electron repulsion? Type yes or no: yes
Give the angular frequency: 0.5
```


### main.cpp
Contains the main program to solve the eigenvalue problems.

### functions.cpp and functions.h
Contains all the necessary functions for main.cpp to work. It's essentially a package deal.

### make_plot.py
A python script that reads data and creates a plot of a wavefunction.


## Results

The [results](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project2/codes/results) are organized as follows:
1. [Buckling beam](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project2/codes/results/BucklingBeam)
2. [Quantum dot with one electron](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project2/codes/results/QM_OneElectron)
3. [Quantum dot with two electron with electron-electron repulsion](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project2/codes/results/QM_TwoElectrons/Repulsion)
4. [Quantum dot without electron-electron repulsion](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project2/codes/results/QM_TwoElectrons/NoRepulsion)
