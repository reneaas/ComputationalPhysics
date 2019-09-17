#include <armadillo>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace arma;


int main(int argc, char* argv[]){
  int n = atoi(argv[1]);
  mat A = mat(n,n);
  A(0,n-1) = 25.0;
  A.print("A = ");
  return 0;
}
