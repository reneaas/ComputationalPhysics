#include <cmath>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "functions.h"

using namespace std;
ofstream ofile;


int main(int nargs, char* args[]){
  int d = atoi(args[1]);        //The spatial dimension of the problem.

  //Solves the 1 + 1 dimensional diffusion equation.
  if(d == 1){
    //Declaration of variables.
    double **v, *t, *x;
    int timesteps, gridpoints, time_index;
    double r, dt, dx, total_time;
    double start_x, end_x;
    string method, outfilename;

    //Command line arguments
    dx = atof(args[2]);
    method = string(args[3]);
    outfilename = string(args[4]);
    time_index = atoi(args[5]);

    //Hardcode variables.
    start_x = 0.;
    end_x = 1.;
    r = 0.5;
    total_time = 0.1;
    dt = r*dx*dx;
    gridpoints = int((end_x - start_x)/dx - 2);
    timesteps = int(total_time/dt - 1);

    //Initiate empty solution matrix.
    v = new double*[timesteps];
    for (int i = 0; i < timesteps; i++){
      v[i] = new double[gridpoints];
    }
    for (int i = 0; i < timesteps; i++){
      for (int j = 0; j < gridpoints; j++){
        v[i][j] = 0.;
      }
    }

    x = new double[gridpoints];
    t = new double[timesteps];

    for (int i = 0; i < gridpoints; i++) x[i] = dx*(i+1);           //Position array.
    for (int i = 0; i < timesteps; i++) t[i] = dt*i;                //Time array
    for (int i = 0; i < gridpoints; i++) v[0][i] = -x[i];           //Initial condition

    if (method == "explicit"){
      Explicit_scheme(v, r, gridpoints, timesteps);
      for (int m = 0; m < timesteps; m++){
        for (int j = 0; j < gridpoints; j++){
          v[m][j] += x[j];
        }
      }


      printf("Writing to file for time-element t[%d]. \n", time_index);
      ofile.open(outfilename);
      ofile << t[time_index] << endl;
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v[time_index][i] << endl;
      }
      ofile.close();

    }

    if (method == "implicit"){
      double *a, *b, *c;
      double *y;

      a = new double[gridpoints];
      b = new double[gridpoints];
      c = new double[gridpoints];
      y = new double[gridpoints];

      for (int m = 0; m < timesteps - 1; m++){
        for (int i = 0; i < gridpoints; i++){
          a[i] = -r;
          b[i] = 1.0 + 2*r;
          c[i] = -r;
          y[i] = v[m][i];
        }
        Forward_substitution(a, b, c, y, gridpoints);
        Back_substitution(v[m+1], b, c, y, gridpoints);
      }
      for (int m = 0; m < timesteps; m++){
        for (int j = 0; j < gridpoints; j++){
          v[m][j] += x[j];
        }
      }


      printf("Writing to file for time-element t[%d]. \n", time_index);
      ofile.open(outfilename);
      ofile << t[time_index] << endl;
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v[time_index][i] << endl;
      }
      ofile.close();
    }


    if (method == "CN"){
      double *a, *b, *c, *y;
      double alpha, beta, gamma;
      a = new double[gridpoints];
      b = new double[gridpoints];
      c = new double[gridpoints];
      y = new double[gridpoints];

      //Hardcode variables
      alpha = r; beta = 1 - 2*r; gamma = r;
      for (int m = 0; m < timesteps - 1; m++){
        for (int i = 0; i < gridpoints; i++){
            a[i] = -r;
            b[i] = 1.0 + 2*r;
            c[i] = -r;
            if (i == 0){
              y[i] =  beta*v[m][i] + gamma*v[m][i+1];
            }
            else if (i == gridpoints-1){
              y[i] =  beta*v[m][i] + alpha*v[m][i-1];
            }
            else{
              y[i] = beta*v[m][i] + gamma*v[m][i+1] + alpha*v[m][i-1];
            }
        }
        //Solving the matrix equation Mv[m+1] = Nv[m] = y
        Forward_substitution(a, b, c, y, gridpoints);
        Back_substitution(v[m+1], b, c, y, gridpoints);
      }



      for (int m = 0; m < timesteps; m++){
        for (int j = 0; j < gridpoints; j++){
          v[m][j] += x[j];
        }
      }

      printf("Writing to file for time-element t[%d]. \n", time_index);
      ofile.open(outfilename);
      ofile << t[time_index] << endl;
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v[time_index][i] << endl;
      }
      ofile.close();
    }
  }
  return 0;
}
