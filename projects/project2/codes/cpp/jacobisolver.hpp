#ifndef JACOBISOLVER_HPP
#define JACOBISOLVER_HPP

#include <armadillo>
#include <iostream>
#include <fstream>

using namespace arma;
using namespace std;

class JacobiSolver {
protected:
    int m_N;
    mat m_A;
    vec m_x, m_eigenvalues;
    double m_max_matrix_val, m_s, m_c, m_t;
    int m_p, m_q;
    ofstream m_ofile;
public:
    JacobiSolver(int N, double a, double b);  //Constructor with V(x) = 0.
    JacobiSolver(int N, double a, double b, double V(double x));  //More general constructor with arbitrary V(x).
    void find_max_off_diagonal();
    void schur2();
    void rotate();
    void compute_eigenvalues(double epsilon);
    void print_eigenvalues();
    void find_ground_state();
};


#endif
