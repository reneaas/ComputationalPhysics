# Project 2: Eigenvalue problems



## Codes
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

3. Quantum dots with two-electrons:

```console
python3 main.py Number_of_gridpoints Max_number_of_iterations qm2
```

* Here you will be prompted to insert additional values for :
1. Repulsion: Whether to include electron repulsion or not.
2. Angular frequency: A double floating point number.

On every run, you'll be prompted with whether you want to compile the code again. If no changes has been made to the code, then type "no". Otherwise type "yes".

### main.cpp
Contains the main program to solve the eigenvalue problems.

### functions.cpp and functions.h
Contains all the necessary functions for main.cpp to work. It's essentially a package deal.
