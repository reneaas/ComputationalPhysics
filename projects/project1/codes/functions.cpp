/*
Here we define the content of the functions declared in the header file.
To include this in your compilation, you can write
1. "c++ -O3 -c project1.cpp functions.cpp"
2. "c++ -O3 -o project1.exe project1.o functions.o"
*/


#include <cmath>      //included some of the function descriptions below used functions from <cmath>.
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;

void LU_thomas(double *a, double *b, double *c, double *d, double *l, double *u, double *y, double *q, double *v, int n){
  //Step 1 and 2: LU-decomposition and forward substitution
  for (int i = 0; i < n; i++){
    if (i == 0){
      d[i] = b[i];
      u[i] = c[i];
      y[i] = q[i];
    }
    else{
      l[i] = a[i-1]/d[i-1];
      d[i] = b[i] - l[i]*u[i-1];
      u[i] = c[i];
      y[i] = q[i] - l[i]*y[i-1];
    }
  }
  //No more use for a, b and c so we deallocate their memory here.
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] l;


  //Step 3: Backward-substitution
  for (int i = n-1; i >= 0; i--){
    if (i == n-1){
      v[i] = y[i]/d[i];
    }
    else{
      v[i] = (y[i] - u[i]*v[i+1])/d[i];
    }
  }

  //Deallocates y and d as their no longer needed.
  delete[] y;
  delete[] d;
  delete[] u;
  return;

}
//LU-decomposition of a tridiagonal matrix.
void LU_decomposition(double* a, double* b, double* c, double* d, double* l, double* u, int n){
  for (int i = 0; i < n; i++){
    if (i == 0){
      d[i] = b[i];
      u[i] = c[i];
    }
    else{
      l[i] = a[i-1]/d[i-1];
      d[i] = b[i] - l[i]*u[i-1];
      u[i] = c[i];
    }
  }
  //No more use for a, b and c so we deallocate their memory here.
  delete[] a;
  delete[] b;
  delete[] c;



  return;
}


//Forward substitution using the LU-decomposition from LU_decomposition.
void Forward_substitutionLU(double* y, double* q, double* l, int n){
  for (int i = 0; i < n; i++){
    if (i == 0){
      y[i] = q[i];
    }
    else{
      y[i] = q[i] - l[i]*y[i-1];
    }
  }

  //q and l has served its purpose and is thus deallocated.
  //delete[] q;
  //delete[] l;
  return;
}

//Back-substitution, works in succession with Forward_substitutionLU.
void Back_substitutionLU(double* v, double* y, double* u, double* d, int n){
  for (int i = n-1; i >= 0; i--){
    if (i == n-1){
      v[i] = y[i]/d[i];
    }
    else{
      v[i] = (y[i] - u[i]*v[i+1])/d[i];
    }
  }

  //Deallocates y and d as their no longer needed.
  delete[] y;
  delete[] d;
  delete[] u;
  return;
}


//Generalized forward substitution of a tridiagonal matrix eq. Ax = b
void Forward_substitution(double* a, double* b, double* c, double* y, int n){
  for (int i = 1; i < n; i++){
    b[i] -= a[i-1]*c[i-1]/b[i-1];
    y[i] -= a[i-1]*y[i-1]/b[i-1];
  }
  delete[] a;
  return;
}

//Generalized back substitution of a tridiagonal matrix eq. Ax = b.
void Back_substitution(double* x, double* b, double* c, double* y, int n){
  for (int i = n-1; i >= 0; i--){
    if (i == n-1){
      x[i] = y[i]/b[i];
    }
    else{
      x[i] = (y[i]-c[i]*x[i+1])/b[i];
    }
  }

  delete[] b;
  delete[] c;
  delete[] y;
  return;
}

//Specialized algorithm to solve the matrix eq Ax = b of a tridiagonal matrix
//where all the elements of the diagonal are equal, and the off-diagonal elements are equal.
/*
void SpecialThomas(double* y, double* x, int n){
  for (int i = 1; i < n; i++){
    y[i] = y[i] + (((double) i)/( (double) (i+1)))*y[i-1];
    //cout << y[i]<< endl;
  }
  for (int i = n-1; i >= 0; i--){
    if (i == n-1){
      x[i] = (double) ((n+1)/n)*y[i];
    }
    else{
      x[i] = (((double)(i+1))/((double)(i+2)))*(y[i] + x[i+1]);
    }
  }
  delete[] y;
  return;
}
*/

void SpecialThomas(double* y, double* x, int n){
  for (int i = 1; i < n; i++){
    y[i] = y[i] + (((double) i)/( (double) (i+1)))*y[i-1];
    //cout << y[i]<< endl;
  }
  for (int i = n-1; i >= 0; i--){
    if (i == n-1){
      x[i] = (((double) n)/((double) (n+1)))*y[i];
    }
    else{
      x[i] = (((double) (i+1))/((double) (i+2)))*(y[i] + x[i+1]);
    }
  }
  delete[] y;
  return;
}

//RHS of the differential eq.
void f(double x, double h, double& vector_element){
   vector_element = 100*exp(-10*x)*h;
   return;
}

//Closed form solution of the differential equation
void closed_form_solution(double x, double& DE_solution){
  DE_solution = 1 - (1 - exp(-10))*x - exp(-10*x);
  return;
}

//Computes the relative error between the closed form solution and the numerical approximation computed in the program project1.cpp.
void compute_errors(double* errors, double* v, double* DE_solution, int n){
  for (int i = 0; i < n; i++){
    errors[i] = log10(abs((v[i] - DE_solution[i])/DE_solution[i]));
  }
  return;
}
