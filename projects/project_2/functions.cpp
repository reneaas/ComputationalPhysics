//Specification of functions.

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <string>

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

string OrthonormalityPreservationTest(mat A, mat S, vec initial_eigval, int n){
  /*
  Units test that check that Orthonormality is preserved by the matrix unitary rotation matrix S.
  */
  string message;
  double norm;
  double tolerance = 1e-4;
  vec eigval;
  mat eigvec;
  vec eigvec1, eigvec2;
  eigvec1 = vec(n);
  eigvec2 = vec(n);
  eig_sym(eigval, eigvec, A);
  for (int i = 0; i < n; i++){
    //eigvec(i) = normli;
  }

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      for (int k = 0; k < n; k++){
        eigvec1(k) = eigvec(i,k);
        eigvec2(k) = eigvec(j,k);
      }
      //norm = dot(eigvec1,eigvec2);
      norm = dot(trans(S*eigvec1), (S*eigvec2));
      if (i == j &&  abs(norm - 1) >= tolerance){
        message = "Ortonormality is not preserved., An eigenvector is not normalized to unity";
        cout << "norm = " << norm << endl;
        return message;
      }
      if ( i != j && abs(norm-0) >= tolerance){
        message = "Orthonormality is not preserved. A pair eigvec aren't orthogonal.";
        cout << "norm = " << norm << endl;
        return message;
      }
    }
  }
  message = "OK";

  return message;
}
