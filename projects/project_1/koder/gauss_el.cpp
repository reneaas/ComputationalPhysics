#include <iostream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char * argv[]){
  mat A =  mat(3,3);
  A.fill(1.0);
  A(0,0) = 1.0;
  A(0,1) = 2.0;
  A(0,2) = -1.0;
  A(1,0) = 2.0;
  A(1,1) = 3.0;
  A(1,2) = -3.0;
  A(2,0) = -1.0;
  A(2,1) = 2.0;
  A(2,2) = 3.0;
  cout << A << endl;
  cout << det(A)<< endl;
  vec b = vec(3);
  b.fill(3.0);
  cout << b << endl;

  for (int m = 0; m < 3; m++){
    for (int j = m+1; j < 3; j++){


      for (int k = m; k < 3; k++){
        A(j,k) = A(j,k) - (A(j,m)*A(k,m))/A(m,m);
          }
      }
  }
  cout << A  << endl;
  cout << b << endl;
  return 0;
}
