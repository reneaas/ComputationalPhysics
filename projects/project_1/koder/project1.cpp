//Project 1 in Computational Physics - FYS3150 - René Ask
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;


double u(double); // Declaration of closed form solution function for project 1.
double ** CreateMatrix(int, int);
void DestroyMatrix(double**, int, int);

int main(int argc, char* argv[]){
  //Using dynamic memory allocation to create matrix A.
  int n = atoi(argv[1]); //First argument in command line is number of points in the n x n matrix.
  double**A = CreateMatrix(n,n);
  for (int i = 0; i < n; i++){
    A[i][i] = 2.0;
  }
  for (int i = 0; i < n-1; i++){
    A[i+1][i] = -1.0;
    A[i][i+1] = -1.0;
  }




  /*
  //General implementation of the tridiagonal matrix using armadillo as framework
  mat A = mat(n,n);
  A.fill(0.0);
  for (int i = 0; i < n; i++){
    A(i,i) = 2.0;
  }
  for (int i = 0; i < n-1; i++){
    A(i+1,i) = -1.0;
    A(i,i+1) = -1.0;
  }
  cout << A << endl;

  //Finding a LU-decomposition using armadillo without the permutation matrix P.
  mat L,U;
  lu(L,U,A);
  cout << "L = " << L << endl;
  cout << "U = " << U << endl;
  */
  return 0;
}


//function specifications
double u(double x){
  // Closed form solution of the differential equations in project 1.
  return 1 - (1-exp(-10))*x - exp(-10*x);
}

double ** CreateMatrix(int n, int m){
  //Allocates memory for a double type matrix and fills it with zeros
  double**mat = new double*[n];
  for (int i = 0; i < n; i++){
    mat[i] = new double[m];
  }
  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      mat[i][j] = 0.0;
    }
  }
  return mat;
}

void DestroyMatrix(double**mat, int n, int m){
  //Deallocated memory of a double type matrix
  for (int i = 0; i < n; i++){
    delete[] mat[i];
  }
  delete[] mat;
}
