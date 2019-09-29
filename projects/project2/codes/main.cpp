#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include "functions.h"
#include <string>
#include "time.h"

using namespace std;
using namespace arma;
ofstream ofile_eigenvalues, ofile_wavefunction, ofile_NumberOfIterations, ofile_Time;

//Main program

int main(int argc, char* argv[]){

  //Declaration of variables:
  int n, N, RowIndex, ColumnIndex, k, l, max_iterations;                      //Integers
  double sinus, cosinus, tangens, tau, tolerance, h, a, d, max_element, rho_max;      //Floating points.
  mat A, S;      //Matrices.
  string message;
  char *outfilename_eigenvalues, *outfilename_wavefunction, *outfilename_NumberOfIterations, *outfilename_Time;
  string problemtype;


  //Specify integers:
  //n = atoi(argv[1]);
  //cin >> n; //Temporary solution to specify the number n from "terminal".
  n = atoi(argv[1]);
  N = n+1;
  max_iterations = atoi(argv[2]);
  outfilename_eigenvalues = argv[3];
  problemtype = string(argv[4]);


  //Specify floats:
  tolerance = 1e-12;
  max_element = 1.0;    //initial value to pass first check in while loop.


  //Create matrices:
  A = mat(n,n);
  S = mat(n,n);
  A.fill(0.0);
  S.fill(0.0);


  //This part fills the matrix A depending on which problem to solve.

  if (problemtype == "BucklingBeam"){
    outfilename_NumberOfIterations = argv[5];
    outfilename_Time = argv[6];
    h = 1.0/((double) N);
    d = 2.0/(h*h);
    a = -1.0/(h*h);
    //Fill up the tridiagonal matrix A:
    for (int i = 0; i < n; i++){
      if (i < n-1){
        A(i,i) = d;
        A(i,i+1) = a;
        A(i+1,i) = a;
      }
      else{
        A(i,i) = d;
      }
    }
  }

  if (problemtype == "QM_OneElectron"){
    rho_max = atof(argv[5]);
    //double rho_max = 8;
    h = rho_max/((double) N);
    d = 2.0/(h*h);
    a = -1.0/(h*h);
    double rho;
    //Now A plays the role of the Hamiltonian matrix
    for (int i = 0; i < n; i++){
      if (i < n-1){
        rho = (i+1)*h;
        A(i,i) = d + rho*rho;
        A(i,i+1) = a;
        A(i+1,i) = a;
      }
      else{
        A(i,i) = d + (i+1)*h;
      }
    }
  }

  if (problemtype == "QM_TwoElectrons"){
    double oscillator_frequency;
    oscillator_frequency = atof(argv[5]);
    outfilename_wavefunction = argv[6];
    string repulsion = string(argv[7]);
    double rho_max = 30;
    h = rho_max/((double) N);
    d = 2.0/(h*h);
    a = -1.0/(h*h);
    double rho;
    /*
    The potential energy is different depending on repulsion is included or not. Therefore
    the Hamiltonian matrix A is computed differently depending on inclusion of repulsion or not.
    */
    if (repulsion == "yes"){
      for (int i = 0; i < n; i++){
        if (i < n-1){
          rho = (i+1)*h;
          A(i,i) = d + rho*rho*oscillator_frequency*oscillator_frequency + 1/rho;
          A(i,i+1) = a;
          A(i+1,i) = a;
        }
        else{
          rho = (i+1)*h;
          A(i,i) = d + rho*rho*oscillator_frequency*oscillator_frequency + 1/rho;
        }
      }
    }
    else{
      for (int i = 0; i < n; i++){
        if (i < n-1){
          rho = (i+1)*h;
          A(i,i) = d + rho*rho*oscillator_frequency*oscillator_frequency;
          A(i,i+1) = a;
          A(i+1,i) = a;
        }
        else{
          rho = (i+1)*h;
          A(i,i) = d + rho*rho*oscillator_frequency*oscillator_frequency;
        }
      }
    }
  }


  //Measuring time for armadillo's eigenvalue function

  time_t start, finish;

  start = clock();

  vec eigval = eig_sym(A);

  finish = clock();

  double time_armadillo = (double) (finish - start)/ (CLOCKS_PER_SEC);
  cout<<"Time used for Armadillo's eigenvalue method:"<<endl;
  cout << "Time used for n = "<< n << ": " << time_armadillo << " s" << endl;



  //Here we find the initial eigenvalues and eigenvectors of the matrix H.
  vec initial_eigenvalues;
  mat initial_eigenvectors;
  eig_sym(initial_eigenvalues, initial_eigenvectors, A);

  int iterations = 0;

  //Measuring time for the Jacobi method

  start = clock();

  //Main algorithm - Jacobi's method
  while (max_element*max_element >= tolerance && iterations < max_iterations){
    Find_MaxElement_and_MaxIndices(RowIndex, ColumnIndex, max_element, A, n);
    Compute_Trigonometric_Functions(RowIndex, ColumnIndex, A, n, tau, tangens, cosinus, sinus);
    k = RowIndex;
    l = ColumnIndex;
    iterations += 1;
    //cout << "iteration = " << iterations << endl;

    //Different version of the algorithm
    double a_kk = A(k,k);
    double a_ll = A(l,l);
    double a_kl = A(k,l);
    A(k,k) = cosinus*cosinus*a_kk - 2.0*cosinus*sinus*a_kl + sinus*sinus*a_ll;
    A(l,l) = sinus*sinus*a_kk + 2.0*cosinus*sinus*a_kl + cosinus*cosinus*a_ll;
    A(l,k) = 0.0;
    A(k,l) = 0.0;

    for (int i = 0; i < n; i++){
      if ( i != k && i != l){
        double a_ik = A(i,k);
        double a_il = A(i,l);
        A(i,k) = cosinus*a_ik - sinus*a_il;
        A(k,i) = A(i,k);
        A(i,l) = cosinus*a_il + sinus*a_ik;
        A(l,i) = A(i,l);
      }
    }
  }

  finish = clock();


  double time_jacobi = (double) (finish - start)/ (CLOCKS_PER_SEC);
  cout<<"Time used for Jacobi's method:"<<endl;
  cout << "Time used for n = "<< n << ": " << time_jacobi << " s" << endl;

  ofile_Time.open(outfilename_Time);
  ofile_Time<<"Jacobi_time"<<" " <<time_jacobi<< " "<<"Armadillo_time" << " "<<time_armadillo<<endl;
  ofile_Time.close();


  S = FillUnitaryMatrix(k, l, n, cosinus, sinus);


  //Unit tests to check if mathematical properties are conserved.
  message = OrthonormalityPreservationTest(A, S, n);
  if (message != "OK"){
    cout << message << endl;
    exit(1);
  }

  //Unit test to check if eigenvalues of matrix A are conserved through the unitary transformation(s).
  message = ConservationOfEigenvalues(A, initial_eigenvalues, n);
  if (message != "OK"){
    cout << message << endl;
    exit(2);
  }


  vec computed_eigenvalues = A.diag();                                          //Extract the computed eigenvalues from the diagonal of the final similar matrix.
  cout << "Completed computations in " << iterations << " iterations" << endl;

  //Here we write the computed eigenvalues to a file.

  computed_eigenvalues = sort(computed_eigenvalues, "ascend");
  ofile_eigenvalues.open(outfilename_eigenvalues);
  ofile_eigenvalues<<rho_max<<endl;
  for (int i = 0; i < n; i++){
    ofile_eigenvalues << setw(20) << setprecision(14) << computed_eigenvalues(i) << endl;
  }
  ofile_eigenvalues.close();


  //If we're solving the two electron problem, we compute the ground state wavefunction and write it to a file.

  if (problemtype == "QM_TwoElectrons"){
    vec ground_state = vec(n);                          //GS wavefunction
    for (int i = 0; i < n; i++){
      ground_state(i) = initial_eigenvectors(0,i);      //Fills the GS wavefunction.
    }

    //Writes the GS wavefunction to file.
    double rho;
    ofile_wavefunction.open(outfilename_wavefunction);
    for (int i = 0; i < n; i++){
      rho = (i+1)*h;
      ofile_wavefunction << ground_state(i)*ground_state(i) << " ";       //It's squared because we want the probability density.
      ofile_wavefunction << rho << endl;
    }
    ofile_wavefunction.close();
  }

  if (problemtype == "BucklingBeam"){
    ofile_NumberOfIterations.open(outfilename_NumberOfIterations);
    ofile_NumberOfIterations << "n" << " " << "iterations" << endl;
    ofile_NumberOfIterations << n << " " << iterations << endl;
    ofile_NumberOfIterations.close();
  }


  return 0;
}
