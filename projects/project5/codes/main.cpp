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
      Explicit_scheme_1D(v, r, gridpoints, timesteps);
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
      double *q;

      a = new double[gridpoints];
      b = new double[gridpoints];
      c = new double[gridpoints];
      q = new double[gridpoints];

      for (int m = 0; m < timesteps - 1; m++){
        for (int i = 0; i < gridpoints; i++){
          a[i] = -r;
          b[i] = 1.0 + 2*r;
          c[i] = -r;
          q[i] = v[m][i];
        }
        Forward_substitution(a, b, c, q, gridpoints);
        Back_substitution(v[m+1], b, c, q, gridpoints);
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
      double *a, *b, *c, *q;
      double alpha, beta, gamma, rho;
      a = new double[gridpoints];
      b = new double[gridpoints];
      c = new double[gridpoints];
      q = new double[gridpoints];
      rho = r/2;                              //In the Crank-Nicolson scheme, r is half of the r used in the two other 1D schemes.

      //Hardcode variables
      alpha = rho; beta = 1 - 2*rho; gamma = rho;
      for (int m = 0; m < timesteps - 1; m++){
        for (int i = 0; i < gridpoints; i++){
            a[i] = -rho;
            b[i] = 1.0 + 2*rho;
            c[i] = -rho;
            if (i == 0){
              q[i] =  beta*v[m][i] + gamma*v[m][i+1];
            }
            else if (i == gridpoints-1){
              q[i] =  beta*v[m][i] + alpha*v[m][i-1];
            }
            else{
              q[i] = beta*v[m][i] + gamma*v[m][i+1] + alpha*v[m][i-1];
            }
        }
        //Solving the matrix equation Mv[m+1] = Nv[m] = y
        Forward_substitution(a, b, c, q, gridpoints);
        Back_substitution(v[m+1], b, c, q, gridpoints);
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


if (d == 2){

  double ***v, *t, *x, *y;
  int timesteps, gridpoints, time_index;
  double r, dt, h, total_time;
  double start_x, end_x;
  string method, outfilename;
  double a,b;

  a = 0.; b = 1;


  //Command line arguments
  h = atof(args[2]);
  outfilename = string(args[3]);
  time_index = atoi(args[4]);

  //Hardcode variables.
  start_x = 0.;
  end_x = 1.;
  r = 0.25;
  total_time = 0.1;
  dt = r*h*h;
  gridpoints = int((end_x - start_x)/h - 2);
  timesteps = int(total_time/dt - 1);

  //Initiate empty solution matrix.
  v = new double**[timesteps];
  for (int i = 0; i < timesteps; i++){
    v[i] = new double*[gridpoints];
    for (int j = 0; j < gridpoints; j++){
      v[i][j] = new double[gridpoints];
    }
  }

  for (int m = 0; m < timesteps; m++){
    for (int i = 0; i < gridpoints; i++){
      for (int j = 0; j < gridpoints; j++){
        v[m][i][j] = 0.;
      }
    }
  }

  x = new double[gridpoints];
  y = new double[gridpoints];
  t = new double[timesteps];

  for (int i = 0; i < gridpoints; i++){
    x[i] = h*(i+1);                                             //Position array x-direction
    y[i] = h*(i+1);                                             //Position array y-direction
  }

  for (int i = 0; i < timesteps; i++) t[i] = dt*i;                //Time array


  //Initial condition
  for (int i = 1; i < gridpoints-1; i++){
    for (int j = 1; j < gridpoints-1; j++){
      v[0][i][j] = (a-b)*y[j] - a;
    }
  }


  Explicit_scheme_2D(v, r, gridpoints, timesteps);




  for (int m = 0; m < timesteps; m++){
    for (int i = 0; i < gridpoints; i++){
      for (int j = 0; j < gridpoints; j++){
        v[m][i][j] += (b-a)*y[j] + a;
      }
    }
  }

  printf("Writing to file for time-element t[%d]. \n", time_index);
  ofile.open(outfilename);
  ofile << t[time_index] << endl;
  for (int i = 0; i < gridpoints; i++){
    for (int j = 0; j < gridpoints; j++){
    //ofile << x[i] << " " << y[j] << " " << v[time_index][i][j] << endl;
    ofile << v[time_index][i][j] << " ";
    }
    ofile << " " << endl;
  }
  ofile.close();


  }

  return 0;
}
