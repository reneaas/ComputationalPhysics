//Project 1 in Computational Physics - FYS3150 - René Ask
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;


double u(double); // Declaration of closed form solution function for project 1.
double f(double); //Declaration of the RHS of the differential equation to be solved.
double ** CreateMatrix(int, int);
void DestroyMatrix(double**, int, int);

int main(int argc, char* argv[]){
  //Using dynamic memory allocation to create matrix A.
  int n = atoi(argv[1]); //First argument in command line is number of points in the n x n matrix.
  double h = atof(argv[2]);

  /*
  double**A = CreateMatrix(n,n);
  for (int i = 0; i < n; i++){
    A[i][i] = 2.0;
  }
  for (int i = 0; i < n-1; i++){
    A[i+1][i] = -1.0;
    A[i][i+1] = -1.0;
  }

  //Gaussian elimination of the generalized tridiagonal matrix:
  for (int m = 0; m < n; m++){
    for (int j = m+1; j < m+3; j++){
      for (int k = m+1; k < m+3; k++){
        A[j][k] = A[j][k] - (A[j][m]*A[m][k])/A[m][m];
      }
    }
  }
  */


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
  cout << "intial A = " << A << endl;
  double* a = new[n];
  double* b = new[n];
  double* c = new[n-1];
  double* d = new[n];
  double* v = new[n]; //The solution to the differential equation - v(x) \approx u(x).
  //Fill the d-array with function values f(x)
  for (int i = 0; i < n; i++){
    double x = ((double)i)*h
    d[i] = f(x)*h*h;
    b[i] = 2.0;
  }
  for (int i = 0; i < n-1; i++){
    c[i] = -1.0;
  }
  for (int i = 1; i < n; i++){
    a[i] = -1.0;
  }

  for (int i = 0; i < n-1; i++){
    if (i = 0){
      c[i] = c[i]/b[i];
      d[i] = d[i]/b[i];
    }
    else{
      c[i] = c[i]/(b[i]-a[i]c[i-1]);
      d[i] = (d[i] - a[i]*d[i-1])/(b[i]-a[i]c[i-1]);
    }
  }

  //Computing solution:
  v[n-1] = d[n-1];
  int i = n-2;
  while (i > 0){
    v[i] = d[i] - c[i]*v[i+1];
  }
  

  cout << A << endl;
  /*
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

double f(double x){
  //RHS of the differential equation
  return 100*exp(-10*x);
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
