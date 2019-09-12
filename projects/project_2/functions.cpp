#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <armadillo>

void compute_squared_elements(mat matrix, mat A, int n){
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (i != j){
        matrix(i,j) = A(i,j)*A(i,j)
      }
      else{
        matrix(i,j) = 0.0;
      }
    }
  }
  return;
}

void find_maxelement_and_maxindices(int &row_index, int &column_index, mat A, int n){
  double max_element = 0.0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; i++){
      if (abs(A(i,j)) > max_element && i != j){
        row_index = i;
        column_index = j;
        max_element = abs(A(i,j));
      }
    }
  }
  return;
}

void compute_trigonometric_functions(int row_index, int column_index, mat A, int n, double &tau, double &tan, double &cos, double &sin){
  tau = (A(row_index, row_index) - A(column_index, column_index))/(2*A(row_index, column_index));

  double t1 = -tau + sqrt(1+tau*tau);
  double t2 = -tau - sqrt(1+tau*tau)
  if (t1 > t2){
    tan = t2;
  }
  else{
    tan = t1;
  }
  cos = 1/sqrt(1+t*t);
  ssin = tan*cos;
  return;
}

void create_similarity_transformation(mat S, int row_index, int column_index){
  S.eye();      //First make S the identity matrix
  S(row_index,row_index) =
}
