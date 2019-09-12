#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

//Declaration of functions.

void fill_RHS(double*, double*, int);

int main(int argc, char* argv[]){
  //Declaration of variables
  int n = atoi(argv[1]); //Number of grid points
  double *x, *y, h;     //Declaration of pointers and stepsize.

  h = 1/((double) n + 1.0);     //Specification of stepsize
  x = new double[n];            //x is specified to be array of length n
  y = new double[n];

  //Fills x with points
  for (int i = 0; i < n; i++){
    x[i] =  ((double) i+1)*h;
  }


  fill_RHS(x, y, n);  //Fills the RHS.

  return 0;
}


void fill_RHS(double *x, double *y, int n){
  for (int i = 0; i < n; i++){
    y[i] = 100*exp(-10*x[i]);
    cout << y[i] << endl;
  }
   return;
}
