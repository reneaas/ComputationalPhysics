#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdio.h>

using namespace std;


double f(double); // Declaration of function

int main(int argc, char * argv[]){
  //Declaration of variables
  int n = 5;
  double stepsize = 1;
  double *derivative1 = new double[n];
  double *derivative2 = new double[n];
  double x = atof(argv[1]);
  double*h  = new double[n];
  for (int i = 0; i < n; i++){
    h[i] = stepsize/(pow(2, (double) i));
    //cout << h[i] << endl;
  }

  for (int i = 0; i < n; i++){
    derivative1[i] = (f(x + h[i])-f(x))/(h[i]);
    cout << derivative1[i] << endl;
  }

  delete h; //Deallocate memory
  return 0;
}


double f(double x){
  return atan(x);
}
