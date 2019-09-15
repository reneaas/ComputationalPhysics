#include <iostream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char * argv[]){
  int n = atoi(argv[1]); //First argument in terminal should to specify the size of the matrix.
  mat A =  mat(n,n);
  for (int i = 0; i < n-1; i++){
    A(i,i) = 2.0;
    A(i,i+1) = -1.0;
    A(i+1,i) = -1.0;
  }
  A(n-1,n-1) = 2.0;
  cout << "A = " << A << endl;
  cout << "det(A) = " << det(A) << endl;
  for (int m = 0; m < n; m++){
    for (int j = m+1; j < n; j++){


      for (int k = m; k < n; k++){
        A(j,k) = A(j,k) - (A(j,m)*A(k,m))/A(m,m);
          }
      }
  }
  cout << A  << endl;
  return 0;
}
