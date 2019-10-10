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

double testfunk(double,double);
double probdens(double,double);

int main(){
  double x1,x2,a,b,mu1,mu2;
  int n;
  n = 10;
  a = 0;
  b = 10;
  double integral = 0;
  double *integrals;
  int N;
  N = 1000000;
  integrals = new double[N];

  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomNumberGenerator(0,1);

  for (int i = 0; i < N; i++){
    cout << "Computing for sample = " << i << endl;
    for (int j = 0; j < n; j++){
      x1 = RandomNumberGenerator(gen);
      mu1 = -log(1-x1);
      for (int k = 0; k < n; k++){
        x2 = RandomNumberGenerator(gen);
        mu2 = -log(1-x2);
        integrals[i] += testfunk(mu1,mu2)/probdens(mu1,mu2);
      }
    }
    integral += integrals[i]/(pow((double) n, 2));
  }
  integral /= (double) N;
  cout << "Computed integral = " << integral << endl;


}

double testfunk(double x1, double x2){
  return x1*x1*x2*x2*exp(-x1-x2);
}

double probdens(double x1, double x2){
  return exp(-x1-x2);
}
