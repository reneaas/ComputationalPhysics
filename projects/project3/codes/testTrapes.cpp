#include <iostream>
#include "lib.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <string>
#include "time.h"
#include <random>
#include <mpi.h>

using namespace std;

double func(double);

int main(){

  int N;
  N = 100;
  double integral,h,x,a,b;
  a = 0;
  b = 1;
  h = 1/N;
  integral = 0;
  integral += (func(a)+func(b))/2;
  //#pragma omp parallel for
  for(int i=0; i< N; i++){
    x = a + i*h;
    integral += func(x);

  }

cout<<integral<<endl;

}

double func(double x){
  return x;
}
