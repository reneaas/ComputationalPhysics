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
    double *v_new, *v_old, *x;
    int gridpoints;
    double r, t, dt, dx, total_time;
    double start_x, end_x;
    string method, outfilename;

    //Command line arguments
    dx = atof(args[2]);
    method = string(args[3]);
    outfilename = string(args[4]);
    total_time = atof(args[5]);
    r = atof(args[6]);

    //Hardcode variables.
    start_x = 0.;
    end_x = 1.;
    dt = r*dx*dx;
    t = 0;


    if (method == "explicit"){
      gridpoints = (int) (end_x - start_x)/dx + 1;
      x = new double[gridpoints];
      for (int i = 0; i < gridpoints; i++) x[i] = dx*(i);           //Position array.

    }
    else{
      gridpoints = (int) ((end_x - start_x)/dx - 1);
      x = new double[gridpoints];
      for (int i = 0; i < gridpoints; i++) x[i] = dx*(i+1);           //Position array.
    }


    //Initiate empty solution vector.
    v_new = new double[gridpoints];
    v_old = new double[gridpoints];

    for (int i = 0; i < gridpoints; i++){
        v_new[i] = 0.;
      }


    for (int i = 0; i < gridpoints; i++) v_old[i] = -x[i];           //Initial condition


    if (method == "explicit"){
      Explicit_scheme_1D(v_new, v_old, r, gridpoints, dt, total_time, t);
      for (int j = 0; j < gridpoints; j++){
        v_new[j] += x[j];
      }


      cout << "Writing to file for t = " << t << endl;
      ofile.open(outfilename);
      ofile << t << endl;
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v_new[i] << endl;
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

      while (t < total_time){
        for (int i = 0; i < gridpoints; i++){
          a[i] = -r;
          b[i] = 1.0 + 2*r;
          c[i] = -r;
          q[i] = v_old[i];
        }
        Forward_substitution(a, b, c, q, gridpoints);
        Back_substitution(v_new, b, c, q, gridpoints);

        for (int k = 0; k < gridpoints; k++) v_old[k] = v_new[k];

        t += dt;
      }
        for (int j = 0; j < gridpoints; j++){
          v_new[j] += x[j];
      }


      cout << "Writing to file for t = " << t << endl;
      ofile.open(outfilename);
      ofile << t << endl;
      ofile << 0.0 << " " <<  0.0 << endl;                                    //writes boundary condition at x = 0: u(0,t) = 0
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v_new[i] << endl;
      }
      ofile << 1.0 << " " << 1.0 << endl;                                    //Writes the boundary condition at x = 1: u(1,t) = 1
      ofile.close();
    }


    if (method == "CN"){
      double *a, *b, *c, *q;
      double alpha, beta, gamma, rho;
      string outfilename2;
      outfilename2 = "contour_data.txt";

      a = new double[gridpoints];
      b = new double[gridpoints];
      c = new double[gridpoints];
      q = new double[gridpoints];
      rho = r/2;                              //In the Crank-Nicolson scheme, r is half of the r used in the two other 1D schemes.

      //Hardcode variables
      alpha = rho; beta = 1 - 2*rho; gamma = rho;
      ofile.open(outfilename2);
      while (t < total_time){
        ofile << t << " ";
        for (int k = 0; k < gridpoints; k++) ofile << v_old[k]+x[k] << " ";
        ofile << endl;
        for (int i = 0; i < gridpoints; i++){
            a[i] = -rho;
            b[i] = 1.0 + 2*rho;
            c[i] = -rho;
            if (i == 0){
              q[i] =  beta*v_old[i] + gamma*v_old[i+1];
            }
            else if (i == gridpoints-1){
              q[i] =  beta*v_old[i] + alpha*v_old[i-1];
            }
            else{
              q[i] = beta*v_old[i] + gamma*v_old[i+1] + alpha*v_old[i-1];
            }
        }
        //Solving the matrix equation Mv[m+1] = Nv[m] = y
        Forward_substitution(a, b, c, q, gridpoints);
        Back_substitution(v_new, b, c, q, gridpoints);


        for (int k = 0; k < gridpoints; k++) v_old[k] = v_new[k];


        t += dt;
      }
      ofile.close();
      for (int j = 0; j < gridpoints; j++){
        v_new[j] += x[j];
      }

      cout << "Writing to file for t = " << t << endl;
      ofile.open(outfilename);
      ofile << t << endl;
      ofile << 0.0 << " " << 0.0 << endl;                                    //writes boundary condition at x = 0: u(0,t) = 0
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v_new[i] << endl;
      }
      ofile << 1.0 << " " << 1.0 << endl;                                    //Writes the boundary condition at x = 1: u(1,t) = 1
      ofile.close();
    }
  }


  if (d == 2){

    double **v_new, **v_old, *x, *y;
    int gridpoints;
    double r, t, dt, h, total_time;
    double start_x, end_x;
    string method, outfilename;
    double n, m;

    //Command line arguments
    h = atof(args[2]);
    outfilename = string(args[3]);
    total_time = atof(args[4]);
    r = atof(args[5]);

    //Hardcode variables.
    start_x = 0.;
    end_x = 1.;
    dt = r*h*h;
    gridpoints = (int) ((end_x - start_x)/h) + 1;
    n = 1;
    m = 1;


    //Initiate empty solution matrix.
    v_new = new double*[gridpoints];
    v_old = new double*[gridpoints];
      for (int j = 0; j < gridpoints; j++){
        v_new[j] = new double[gridpoints];
        v_old[j] = new double[gridpoints];
      }


    for (int i = 0; i < gridpoints; i++){
      for (int j = 0; j < gridpoints; j++){
        v_old[i][j] = 0.;
        v_new[i][j] = 0.;
      }
    }


    x = new double[gridpoints];
    y = new double[gridpoints];

    for (int i = 0; i < gridpoints; i++){
      x[i] = (double) h*i;                                             //Position array x-direction
      y[i] = (double) h*i;                                             //Position array y-direction
    }

    //initial condition
    for (int i = 1; i < gridpoints - 1 ; i++){
      for (int j = 1; j < gridpoints -1; j++){
        //v_old[i][j] = exp(-alpha*(fabs(x[i]-0.5) + fabs(y[j]-0.5)));
        v_old[i][j] = sin(n*M_PI*x[i])*sin(m*M_PI*y[j]);
      }
    }

    Explicit_scheme_2D(v_new, v_old, r, gridpoints, dt, total_time, t);


    cout << "Writing to file for t = " << t << endl;
    ofile.open(outfilename);
    ofile << t << endl;
    for (int i = 0; i < gridpoints; i++){
      for (int j = 0; j < gridpoints; j++){
      ofile << v_new[i][j] << " ";
      }
      ofile << " " << endl;
    }
    ofile.close();

    }

}
