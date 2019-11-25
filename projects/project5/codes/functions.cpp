#include <cmath>      //included some of the function descriptions below used functions from <cmath>.
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;

//Algorithm for forward Euler method
void Explicit_scheme(double **v , double r, int gridpoints, int timesteps){

  for (int m = 0; m < timesteps-1; m++){
    for (int j = 0; j < gridpoints; j++){
      if (j == 0){
        v[m+1][j] = (1-2*r)*v[m][j] + r*v[m][j+1];
      }
      else if (j == gridpoints - 1){
        v[m+1][j] = (1-2*r)*v[m][j] + r*v[m][j-1];
      }
      else{
        v[m+1][j] = (1-2*r)*v[m][j] + r*(v[m][j+1] + v[m][j-1]);
      }

    }
  }
}

//Generalized forward substitution of a tridiagonal matrix eq. Ax = b
void Forward_substitution(double* a, double* b, double* c, double* y, int n){
  for (int i = 1; i < n; i++){
    b[i] -= a[i-1]*c[i-1]/b[i-1];
    y[i] -= a[i-1]*y[i-1]/b[i-1];
  }
  return;
}

//Generalized backward substitution of a tridiagonal matrix eq. Ax = b.
void Back_substitution(double* x, double* b, double* c, double* y, int n){
  for (int i = n-1; i >= 0; i--){
    if (i == n-1){
      x[i] = y[i]/b[i];
    }
    else{
      x[i] = (y[i]-c[i]*x[i+1])/b[i];
    }
  }
  return;
}
