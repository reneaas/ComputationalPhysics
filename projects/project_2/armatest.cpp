#include <armadillo>
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace arma;

int main(int argc, char* argv[]){
  int n = atoi(argv[1]);
  mat A = mat(n,n);
  for (int i = 0; i < n; i++){
    if (i < n-1){
      A(i,i) = 2.0;
      A(i,i+1) = -1.0;
      A(i+1,i) = -1.0;
    }
    else{
      A(i,i) = 2.0;
    }
  }
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, A);

  eigvec.print("eigenvectors = ");
  eigval.print("eigenvalues = ");
  cout << "inner product:" << eigvec(0:2) << endl;
}
