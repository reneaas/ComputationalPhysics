#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include "functions.h"


using namespace std;
using namespace arma;

ofstream ofile_eigenvalues;                       //Declaration of Global outputfile to write the computes eigenvalues to.

int main(int argc, char* argv[]){
  //Declaration of variables.
  int n = atoi(argv[1]);                          //Size of matrix in the program
  char* outfilename_eigenvalues;                  //Declaration of outputfile name for the eigenvalues.
  double tolerance = 1e-6;                        //Tolerance for which the while loop will stop once the proper condition is reached.
  double h, a, d;
  h = 1/((double) n);
  d = 2/(h*h);
  a = -1/(h*h);
  int row_index, column_index;
  double tau, tan, cos, sin;
  double max_element = 1.0;

  mat A = mat(n,n);
  mat S = mat(n,n);
  mat S_transpose = mat(n,n);
  mat B = mat(n,n);

  //Main Algorithm
  while (abs(max_element*max_element) > tolerance){
    compute_squared_elements(A_squared_elements, A, n);
    find_maxelement_and_maxindices(row_index, column_index, A, n);
    compute_trigonometric_functions(row_index, column_index, A, n, tau, tan, cos, sin);
    max_element = abs(A(row_index,column_index));

  }
  return 0;
}
