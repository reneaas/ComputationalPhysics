#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <cstdlib>

using namespace std;
using namespace arma;

void CreateMatrix(int m, int n){
  mat = new double*[n];
  for (int i = 0; i < n; i++){
    mat[i] = new double[n];
  }
  return;
}

void DestroyMatrix(double **mat, int n){
  for (int i = 0; i < n; i++){
    delete[] mat[i];
  }
  delete[] mat;
  return;
}

void MatrixMultiplication(double **C, double **A, double **B, int n){
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      for (int k = 0; k < n; k++){
        C[i][j] += A[i][k]*B[k][j];
      }
    }
  }
  return;
}

void FillUnitaryMatrix(double **S, double **S_transpose, int k, int l, int n, double c, double t){
  //First make the matrix the identity matrix.
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i == j){
        S[i][j] = 1.0;
        S_transpose[i][j] = 1.0;
      }
      else{
        S[i][j] = 0.0;
        S_transpose[i][j] = 0.0;
      }
    }
  }
  //Fill in the rest of the elements to make it into a rotation matrix.
  S[k][k] = c;
  S[l][l] = c;
  S[k][l] = -s;
  S[l][k] = s;

  //Fill the transposed unitary matrix.
  S_transpose[k][k] = c;
  S_transpose[l][l] = c;
  S_transpose[l][k] = -s;
  S_transpose[k][l] = s;
  return;
}

void compute_squared_elements(double **matrix, double **A, int n){
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i != j){
        matrix[i,j] = A[i,j]*A[i,j]
      }
      else{
        matrix[i,j] = 0.0;
      }
    }
  }
  return;
}

void find_maxelement_and_maxindices(int &row_index, int &column_index, double &max_element, double **A, int n){
  max_element = 0.0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; i++){
      if (abs(A[i,j]) > max_element && i != j){
        row_index = i;
        column_index = j;
        max_element = abs(A[i,j]);
      }
    }
  }
  return;
}

void compute_trigonometric_functions(int row_index, int column_index, double **A, int n, double &tau, double &tan, double &cos, double &sin){
  tau = (A[row_index, row_index] - A[column_index, column_index])/(2*A[row_index, column_index]);

  double t1 = -tau + sqrt(1+tau*tau);
  double t2 = -tau - sqrt(1+tau*tau)
  if (t1 > t2){
    tan = t2;
  }
  else{
    tan = t1;
  }
  cos = 1/sqrt(1+tan*tan);
  sin = tan*cos;
  return;
}
