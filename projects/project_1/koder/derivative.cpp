#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;

double f(double);   //Declaration

int main(int argc, char* argv[]){
  double derivative;
  double x = atof(argv[1]);
  double h = 0.5;
  derivative = (f(x+h) - f(x))/h;
  cout << derivative < endl;
  return 0;
}


double f(double x){
  return exp(-x);
}
