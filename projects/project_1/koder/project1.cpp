#include <iostream>
#include <cmath>
#include <cstdlib>
#include <armadillo>

using namespace std, arma;


double f(double);

int main(){


  return 0;
}


double f(double x){
  return 1 - (1-exp(-10))*x - exp(-10*x);
}
