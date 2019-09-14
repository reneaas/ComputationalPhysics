//Specification of functions.

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;


void Find_MaxElement_and_MaxIndices(int &RowIndex, int &ColumnIndex, double &max_element, mat A, int n){
  max_element = 0.0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (abs(A(i,j)) > max_element && i != j){
        RowIndex = i;
        ColumnIndex = j;
        max_element = abs(A(i,j));
      }
    }
  }
}


void Compute_Trigonometric_Functions(int row_index, int column_index, mat A, int n, double &tau, double &tangens, double &cosinus, double &sinus){
  double k = row_index;
  double l = column_index;
  if (A(k,l) != 0.0) {
    tau =  ( A(l,l) - A(k,k) ) / (2.0 * A(k,l));
    if (tau >= 0) {
      tangens = 1.0/(tau + sqrt(1.0 + tau * tau));
    } else {
      tangens = -1.0/(-tau + sqrt(1.0 + tau * tau));
    }
    cosinus = 1.0 / sqrt(1.0 + tangens * tangens);
    sinus = tangens*cosinus;
  }
  else{
    cosinus = 1.0;
    sinus = 0.0;
  }
}

mat FillUnitaryMatrix(int k, int l, int n, double cosinus, double sinus){
  mat S = mat(n,n);
  S.eye();
  //Fill in the rest of the elements to make it into a rotation matrix.
  S(k,k) = cosinus;
  S(l,l) = cosinus;
  S(k,l) = -sinus;
  S(l,k) = sinus;
  return S;
}

string OrthonormalityPreservationTest(mat A, mat S, vec initial_eigenvalues, int n){
  string message;
  double norm;
  double tolerance = 1e-4;
  vec eigenvalues;
  mat eigenvectors;
  eig_sym(eigenvalues, eigenvectors, A);
  eigenvectors = normalise(eigenvectors);
  for (int i = 0; i < n; i++){
    //eigenvectors(i) = normli;
  }

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      norm = dot(trans(S*eigenvectors(i)), S*eigenvectors(i));
      if (i == j &&  abs(norm - 1) >= tolerance){
        message = "Ortonormality is not preserved., An eigenvector is not normalized to unity";
        return message;
      }
      if ( i != j && abs(norm-0) >= tolerance){
        message = "Orthonormality is not preserved. A pair eigenvectors aren't orthogonal.";
        return message;
      }
    }
  }
  message = "OK";

  return message;
}
