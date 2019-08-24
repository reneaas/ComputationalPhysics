#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

//Declaration of functions
double** CreateMatrix(int, int);
void DestroyMatrix(double**, int, int);


//Main program
int main(argc, char* argv[]){
  double ** A;  //Double pointer for a matrix. A
  int n = atoi(argv[1]);
  A = CreateMatrix(n,n);
  double*a = new[n];
  double*b = new[n];
  double*c = new[n];
  double*d = new[n];
  double*l = new[n];

  //Create the LU-decomposition of A.
  for (int i = 0; i < n; i++){
    if (i == 1){
      d[0] = a[0];
      u[0] = c[0];
    }

    if (0 < i < n){
      l[i] = b[i]/d[i-1];
      d[i] = a[i] - l[i]*u[i-1];
      u[i] = c[i];
    }

    if (i == n){
      l[n] = b[n]/d[n-1];
      d[n] = a[n] - l[n]*u[n-1];
    }
  }


}



//Specification of functions declared in before main()
double** CreateMatrix(int n, int m){
  //Creates a (n x m)-matrix with dynamic allocation
  double**mat = new*[n];
  for (int i = 0; i < n; i++){
    mat[i] = new[m];
  }

  for (int i = 0; i < n; i++){
    for (int j = 0; j < m; j++){
      mat[i][j] = 0.0;
    }
  }
  return mat;
}

void DestroyMatrix(double** mat, int n, int m){
  //Deletes matrix, that is, it deallocates memory.
  for (int i = 0; i < n; i++){
    delete[] mat[i];
  }
  delete[] mat;
}
