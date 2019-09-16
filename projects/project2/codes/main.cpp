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


//Main program

int main(int argc, char* argv[]){

  //Declaration of variables:
  int n, RowIndex, ColumnIndex, k, l, max_iterations;         //Integers
  double sinus, cosinus, tangens, tau, tolerance, h, a, d, max_element;      //Floating points.
  mat A, B, S;      //Matrices.
  string message;

  //Specify integers:
  //n = atoi(argv[1]);
  //cin >> n; //Temporary solution to specify the number n from "terminal".
  n = atoi(argv[1]);
  max_iterations = atoi(argv[2]);

  //Specify floats:
  h = 1.0/((double) n);
  d = 2.0/(h*h);
  a = -1.0/(h*h);
  tolerance = 1e-8;
  max_element = 1.0;    //initial value to pass first check in while loop.


  //Create matrices:
  A = mat(n,n);
  S = mat(n,n);
  B = mat(n,n);
  B.fill(0.0);
  A.fill(0.0);
  S.fill(0.0);
  //----------------Test zone ahead!--------------------

  //Fill matrices with random values to test function:
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
  //A.print(" Initial matrix A  = ");

  //Here we find the initial eigenvalues and eigenvectors of the matrix A.
  vec initial_eigenvalues;
  mat initial_eigenvectors;
  eig_sym(initial_eigenvalues, initial_eigenvectors, A);
  initial_eigenvalues.print("Eigenvalues of A = ");

  int iterations = 0;

  //Main algorithm
  while (max_element*max_element > tolerance && iterations < max_iterations){
    Find_MaxElement_and_MaxIndices(RowIndex, ColumnIndex, max_element, A, n);
    Compute_Trigonometric_Functions(RowIndex, ColumnIndex, A, n, tau, tangens, cosinus, sinus);
    k = RowIndex;
    l = ColumnIndex;
    S = FillUnitaryMatrix(k, l, n, cosinus, sinus);
    message = OrthonormalityPreservationTest(A, S, n);                          //Units test to check if orthonormality is preserved.
    if (message != "OK"){
      cout << message << endl;
      exit(1);
    }
    A = trans(S)*A*S;                                                           //Computes the similarity matrix, essentially setting A to be the new similar matrix.
    message = ConservationOfEigenvalues(A, initial_eigenvalues, n);             //Units test to check if eigenvalues of matrix A are conserved through the unitary transformation.
    if (message != "OK"){
      cout << message << endl;
      exit(2);
    }
    iterations += 1;
  }
  vec computed_eigenvalues = A.diag();                                          //Extract the computed eigenvalues from the diagonal of the final similar matrix.
  cout << "Results after " << iterations << " iterations" << endl;
  //A.print(" A = ");
  sort(computed_eigenvalues, "ascend").print("computed_eigenvalues = ");                        //Print the computed eigenvalues.
  return 0;
}
