#include "jacobisolver.hpp"

double V(double x);

int main(int argc, char const *argv[]) {
    int N = 250;
    double a = 0.;
    double b = 5.;
    
    JacobiSolver my_solver(N, a, b, V);
    my_solver.compute_eigenvalues(1e-10);
    my_solver.print_eigenvalues();
    return 0;
}


double V(double x){
    double omega = 0.25;
    return omega*omega*x*x + 1./x;
}
