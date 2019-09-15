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
  cout << "Initial matrix A" << endl;
  A.print(" A  = ");

  //Find "analytical eigenvalues and eigenvectors"
  vec initial_eigenvalues;
  mat initial_eigenvectors;
  eig_sym(initial_eigenvalues, initial_eigenvectors, A);
  //normalise(initial_eigenvectors);
  //normalise(initial_eigenvectors);
  initial_eigenvectors.print("Initial eigenvectors =");

  initial_eigenvectors.print("Initial normalised eigenvectors = ");
  //double length = dot(trans(initial_eigenvectors(0)),initial_eigenvectors(0));
  initial_eigenvalues.print("Eigenvalues = ");

  int iterations = 0;
  while (max_element*max_element > tolerance && iterations < max_iterations){
    Find_MaxElement_and_MaxIndices(RowIndex, ColumnIndex, max_element, A, n);
    Compute_Trigonometric_Functions(RowIndex, ColumnIndex, A, n, tau, tangens, cosinus, sinus);
    k = RowIndex;
    l = ColumnIndex;
    S = FillUnitaryMatrix(k, l, n, cosinus, sinus);

    //This part needs to be fixed before it's included in the program.
    //Check if orthonormality is preserved.
    message = OrthonormalityPreservationTest(A, S, n);
    if (message != "OK"){
      cout << message << endl;
      exit(1);
    }

    message = ConservationOfEigenvalues(A, initial_eigenvalues, n);
    if (message != "OK"){
      cout << message << endl;
      exit(2);
    }
    //Then compute the similary matrix and implement the whole charade in a while loop.
    A = trans(S)*A*S;
    //A = B;

    iterations += 1;
  }
  vec computed_eigenvalues = A.diag();
  cout << "Results after " << iterations << " iterations" << endl;
  A.print(" A = ");
  computed_eigenvalues.print("computed_eigenvalues = ");
  return 0;
}
