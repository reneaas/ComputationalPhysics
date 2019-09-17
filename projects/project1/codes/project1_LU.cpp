#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "functions.h"
#include <armadillo>

ofstream ofile_solution;

int main(int argc, char* argv[]){
  int n = atoi(argv[1]);
  char *outfilename_solution, *outfile = argv[2];
  double h = 1/((double) n + 1);
  double hh = h*h;
  mat A = mat(n,n);
  mat L, U;
  vec v, y;
  vec q = vec(n);

  for (int i = 0; i < n; i++){
    if (i < n-1){
    A(i,i) = 2.0;
    A(i,i+1) = -1.0;
    A(i+1,i) = -1.0;
    }
    else{
      A(i,i) = 2.0;
    }
  }

  //Fill the RHS
  for (int i = 0; i < n; i++){
    x = (double) (i+1)*h
    f(x, hh, q(i));
  }

  //LU-decomposition, A = LU.
  lu(L,U,A);


  y = solve(L,q);     //First solve Ly = q
  v = solve(U,y);     //Then solve Uv = y

  ofile.open(outfilename_solution);
  for (int i = 0; i < n; i++){
    ofile << v(i) << endl;
  }
  ofile.close();
}
