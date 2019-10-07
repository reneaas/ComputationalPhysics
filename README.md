# FYS3150 - Computational Physics

This repository is a collaboration between
Kaspara Skovli Gåsvær, Maria Linea Horgen and René Ask. This
contains all our work in the course [FYS3150 - Computational Physics](https://www.uio.no/studier/emner/matnat/fys/FYS3150/)
at UiO.

## Projects
All coursework is found under [projects](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/).

### Project 1: 1D-Possions equation
In [project 1](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project1), the 1-dimensional Poisson equation is solved. The differential equation first recast into a linear tridiagonal matrix equation where the particular matrix we deal with
is known as a Toeplitz matrix. It's solved through three different algorithms: a generalized and a specialized algorithm, and one that
combines the generalized one with LU-decomposition.

Here are the [codes](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project1/codes) and the
[report](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project1/report/project1.pdf)

### Project 2: Eigenvalue problems
In [project 2](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project2) we solve second order differential equations (DEs) by recasting them into matrix equations and solving eigenvalue problems. These are solved by a direct method known as Jacobi's method.

Here are the [codes](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project2/codes) and the
[report](https://github.com/reneaas/ComputationalPhysics/blob/master/projects/project2/report/Project_2_report.pdf)


### Project 3: Integration methods
In [project 3](https://github.com/reneaas/ComputationalPhysics/tree/master/projects/project3), we solve an integral that is used in quantum mechanics to approximate the ground state
of two-electron systems such as the Helium atom. We derive an analytical expression for the integral before we proceed to solve this using two forms of Gaussian quadrature, Gauss-Legendre and Gauss-Laguerre quadrature. We proceed then to solve the same integral using Monte Carlo integration, both with and without importance sampling.
