#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(){
  mat A = randu<mat>(5,5);
  mat B = randu<mat>(5,5);

  cout << A*B << endl;

  mat C;
  C.set_size(3,3);
  mat A_transpose = trans(A);
  A.print("A:");
  A_transpose.print("A^T: ");
  return 0;
}
