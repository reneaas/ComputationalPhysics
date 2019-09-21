#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include "functions.h"
#include <string>

using namespace std;
using namespace arma;
ofstream ofile_eigenvalues, ofile_wavefunction;

//Main program

int main(int argc, char* argv[]){

  //Declaration of variables:
  int n, N, RowIndex, ColumnIndex, k, l, max_iterations;         //Integers
  double sinus, cosinus, tangens, tau, tolerance, h, a, d, max_element;      //Floating points.
  mat A, B, S, Hamiltonian;      //Matrices.
  string message;
  char *outfilename_eigenvalues, *outfilename_wavefunction;
  string problemtype;


  //Specify integers:
  //n = atoi(argv[1]);
  //cin >> n; //Temporary solution to specify the number n from "terminal".
  N = atoi(argv[1]);
  n = N;
  max_iterations = atoi(argv[2]);
  outfilename_eigenvalues = argv[3];
  problemtype = string(argv[4]);


  //Specify floats:
  tolerance = 1e-8;
  max_element = 1.0;    //initial value to pass first check in while loop.


  //Create matrices:
  A = mat(n,n);
  S = mat(n,n);
  B = mat(n,n);
  B.fill(0.0);
  A.fill(0.0);
  S.fill(0.0);

  if (problemtype == "BucklingBeam"){
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
    double rho_max = 30;
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
    //cout << "Specify oscillator frequency" << endl;
    //cin >> oscillator_frequency;
    h = rho_max/((double) N);
    d = 2.0/(h*h);
    a = -1.0/(h*h);
    double rho;
    //Now A plays the role of the Hamiltonian matrix
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



  //A.print(" Initial matrix A  = ");

  //Here we find the initial eigenvalues and eigenvectors of the matrix A.
  vec initial_eigenvalues;
  mat initial_eigenvectors;
  eig_sym(initial_eigenvalues, initial_eigenvectors, A);

  int iterations = 0;

  //Main algorithm
  while (max_element*max_element >= tolerance && iterations < max_iterations){
    Find_MaxElement_and_MaxIndices(RowIndex, ColumnIndex, max_element, A, n);
    Compute_Trigonometric_Functions(RowIndex, ColumnIndex, A, n, tau, tangens, cosinus, sinus);
    k = RowIndex;
    l = ColumnIndex;

    S = FillUnitaryMatrix(k, l, n, cosinus, sinus);
    /*
    message = OrthonormalityPreservationTest(A, S, n);                          //Units test to check if orthonormality is preserved.
    if (message != "OK"){
      cout << message << endl;
      exit(1);
    }
    */
    A = trans(S)*A*S;
    /*
    message = ConservationOfEigenvalues(A, initial_eigenvalues, n);             //Units test to check if eigenvalues of matrix A are conserved through the unitary transformation.
    if (message != "OK"){
      cout << message << endl;
      exit(2);
    }
    */
    cout << "iteration = " << iterations << endl;
    iterations += 1;
  }

  vec computed_eigenvalues = A.diag();                                          //Extract the computed eigenvalues from the diagonal of the final similar matrix.
  cout << "Completed computations in " << iterations << " iterations" << endl;

  //Write the computed eigenvalues to file.

  vec difference_in_eigenvalues = vec(n);
  computed_eigenvalues = sort(computed_eigenvalues, "ascend");
  ofile_eigenvalues.open(outfilename_eigenvalues);
  ofile_eigenvalues << "Analytical eigenvalues:" << endl;
  for (int i = 0; i < n; i++){
    ofile_eigenvalues << initial_eigenvalues(i) << endl;
  }

  ofile_eigenvalues << "Computed eigenvalues:" << endl;
  for (int i = 0; i < n; i++){
    ofile_eigenvalues << computed_eigenvalues(i) << endl;
  }
  ofile_eigenvalues.close();



  if (problemtype == "QM_TwoElectrons"){
    vec ground_state = vec(n);
    for (int i = 0; i < n; i++){
      ground_state(i) = initial_eigenvectors(0,i);
    }

    double rho;
    ofile_wavefunction.open(outfilename_wavefunction);
    for (int i = 0; i < n; i++){
      rho = (i+1)*h;
      ofile_wavefunction << ground_state(i)*ground_state(i) << " ";
      ofile_wavefunction << rho << endl;
    }
    ofile_wavefunction.close();
  }
  return 0;
}
