#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "functions.h"

using namespace std;
ofstream ofile_solution, ofile_errors, ofile_time;


int main(int argc, char* argv[]){
  //Declaration of variables.
  int n = atoi(argv[1]);        //Number of grid points
  char *outfilename_solution, *outfilename_errors, *outfilename_time; //Declaration of filenames
  string which_algorithm;
  outfilename_solution = argv[2];        //The actual filename of the file we'll write v(x) to.
  outfilename_errors = argv[3];          //Filename we'll write the computed relative errors to.
  outfilename_time = argv[4];
  which_algorithm = string(argv[5]);
  double *a, *b, *c, *d, *l, *u, *q, *v, *y, *DE_solution, *errors;    //Declaration of pointers to represent vectors.
  double h;                                     //Declaration of stepsize defined by the start_point, end_point and number of grid points n.
  double start_point = 0.0;
  double end_point = 1.0;
  h = (end_point-start_point)/( (double) n + 1.0);    //Defines the stepsize.


  //Specification of the vectors used in the calculations using dynamic memory allocation.
  a = new double[n];
  b = new double[n];
  c = new double[n];
  d = new double[n];
  l = new double[n];
  u = new double[n];
  q = new double[n];
  v = new double[n];
  y = new double[n];
  DE_solution = new double[n];
  errors = new double[n];

  //Filling the arrays with the necessary values:
  double h_squared = h*h; //Multiplication factor that reduces number of flops in the for-loop.
  for (int i = 0; i < n; i++){
    a[i] = -1.0;
    b[i] = 2.0;
    c[i] = -1.0;
    double x = ((double) i + 1.0)*h;
    f(x, h_squared, q[i]); //Call by reference to speed up execution.
    closed_form_solution(x, DE_solution[i]);
  }


  string general_algorithm = "general_algorithm";
  string special_algorithm = "special_algorithm";
  string LU_Thomas = "LU";

  clock_t start, finish;    //Declaration of variables to time the main algorithm
  start = clock();          //Starts the clock.

  //Main algorithm:
  if (which_algorithm == general_algorithm){
    Forward_substitution(a, b, c, q, n);    //Step 1: Forward substitution
    Back_substitution(v, b, c, q, n);       //Step 2: Back-substitution
    cout << "Running general algorithm" << endl;
  }
  if (which_algorithm == special_algorithm){
    SpecialThomas(q, v, n);
    cout << "Running special algorithm" << endl;
  }
  if (which_algorithm == LU_Thomas){
    cout << "Running LU-Thomas" << endl;
    //LU_decomposition(a,b,c,d,l,u,n);      //Step 1: LU-decomposition - A = LU
    //Forward_substitutionLU(y,q,l,n);      //Step 2: Forward-substitution - solves Ly = q
    //Back_substitutionLU(v,y,u,d,n);       //Step 3: Back-substitution - solves Uv = y
    LU_thomas(a, b, c, d, l, u, y, q, v, n);
  }


  //Compute the time interval the main algorithm took to complete.
  finish = clock();     //Stops the clock
  double timeused = (double) (finish-start)/(CLOCKS_PER_SEC);   //Computes the time elapsed for the main algorithm to finish
  cout << "Total time = " << timeused << " s" << endl;          //Prints the time elapsed to screen.
  //Writes the time elapsed to a file
  ofile_time.open(outfilename_time);
  ofile_time << timeused << endl;
  ofile_time.close();


  //Compute errors
  compute_errors(errors, v, DE_solution, n);

  //Write computed solution and computed relative errors to their own files.
  ofile_solution.open(outfilename_solution);
  ofile_errors.open(outfilename_errors);
  for (int i = 0; i < n; i++){
    ofile_solution << v[i] << endl;
    ofile_errors << errors[i] << endl;
  }
  ofile_solution.close();
  ofile_errors.close();
  delete[] v;
  delete[] errors;


  //Write the time elapsed to a file
  ofile_time.open(outfilename_time);
  ofile_time << timeused << endl;
  ofile_time.close();

  return 0;
}
