#include "jacobisolver.hpp"

JacobiSolver::JacobiSolver(int N, double a, double b){
    m_N = N;
    double h = (b-a)*(1./(N+1));
    double hh_inv = 1./(h*h);
    m_x = linspace(a + h , b - h, m_N);
    double diag_elem = 2./(h*h);

    //Fill matrix
    m_A = mat(m_N, m_N).fill(0.);
    for (int i = 0; i < m_N-1; i++){
        m_A(i,i) = diag_elem;
        m_A(i, i+1) = -hh_inv;
        m_A(i+1, i) = -hh_inv;
    }
    m_A(m_N-1, m_N-1) = diag_elem;
    
}

JacobiSolver::JacobiSolver(int N, double a, double b, double V(double x))
{
    m_N = N;
    double h = (b-a)*(1./(N+1));
    double hh_inv = 1./(h*h);
    m_x = linspace(a + h , b - h, m_N);
    double diag_elem = 2./(h*h);

    //Fill matrix
    m_A = mat(m_N, m_N).fill(0.);
    for (int i = 0; i < m_N-1; i++){
        m_A(i,i) = diag_elem + V(m_x(i));
        m_A(i, i+1) = -hh_inv;
        m_A(i+1, i) = -hh_inv;
    }
    m_A(m_N-1, m_N-1) = diag_elem + V(m_x(m_N-1));
}

void JacobiSolver::find_max_off_diagonal() {

    //Search for max value at the upper off-diagonal of A
    m_max_matrix_val = 0.;
    for (int i = 0; i < m_N; i++){
        for (int j = i+1; j < m_N; j++){
            if (abs(m_A(i,j)) > m_max_matrix_val){
                m_max_matrix_val = abs(m_A(i,j));
                m_p = i;
                m_q = j;
            }
        }
    }
}

void JacobiSolver::schur2() {
    double tau;
    if (m_A(m_p, m_q) != 0){
        tau = ( (m_A(m_q, m_q) - m_A(m_p, m_p))/(2.*m_A(m_p, m_q)) );
        if (tau >= 0){
            m_t = 1./(tau + sqrt(1+tau*tau));
        }
        else{
            m_t = 1./(tau - sqrt(1+tau*tau));
        }
        m_c = 1/sqrt(1. + m_t*m_t);
        m_s = m_t*m_c;
    }
    else{
        m_c = 1.;
        m_s = 0.;
    }
}

void JacobiSolver::rotate(){
    //Perform rotation transform
    double a_pp = m_A(m_p, m_p);
    double a_qq = m_A(m_q, m_q);
    double a_pq = m_A(m_p, m_q);

    m_A(m_p, m_p) = m_c*m_c*a_pp - 2.*m_c*m_s*a_pq + m_s*m_s*a_qq;
    m_A(m_q, m_q) = m_s*m_s*a_pp + 2.*m_c*m_s*a_pq + m_c*m_c*a_qq;
    m_A(m_p, m_q) = 0.;
    m_A(m_q, m_p) = 0.;


    for (int i = 0; i < m_N; i++){
        if (i != m_p && i != m_q){
            double a_ip = m_A(i, m_p);
            double a_iq = m_A(i, m_q);
            m_A(i, m_p) = m_c*a_ip - m_s*a_iq;
            m_A(m_p, i) = m_A(i, m_p);
            m_A(i, m_q) = m_c*a_iq + m_s*a_ip;
            m_A(m_q, i) = m_A(i, m_q);
        }
    }
}

void JacobiSolver::compute_eigenvalues(double epsilon) {
    m_max_matrix_val = 2*epsilon; //Just there to enter the while loop
    int counter = 0;
    while (m_max_matrix_val >= epsilon){
        find_max_off_diagonal();
        schur2();
        rotate();
        find_max_off_diagonal();
        counter++;
    }
    m_eigenvalues = diagvec(m_A); //Extract diagonal elements
    m_eigenvalues = sort(m_eigenvalues, "ascend");
}

void JacobiSolver::print_eigenvalues(){
    cout << "Eigenvalues = " << endl;
    for (int i = 0; i < 10; i++){
        cout << m_eigenvalues(i) << endl;
    }
}

void JacobiSolver::find_ground_state(){
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, m_A);
    vec ground_state = vec(m_N);
    double ground_state_energy = eigval(0);
    for (int i = 0; i < m_N; i++){
        ground_state(i) = eigvec(i, 0);
    }

    string filename = "buckling_beam.txt";
    m_ofile.open(filename);
    for (int i = 0; i < m_N; i++){
        m_ofile << m_x(i) << " " << ground_state(i) << endl;
    }
    m_ofile.close();
}
