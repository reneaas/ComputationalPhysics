#include <cmath>
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

  int d = atoi(args[1]);

  if(d == 1){

    double **v, *t, *x;
    int timesteps, gridpoints;
    double r, dt, dx, total_time;
    double start_x, end_x;
    string method, outfilename;

    dx = atof(args[2]);
    method = string(args[3]);
    outfilename = string(args[4]);

    start_x = 0.;
    end_x = 1.;
    r = 0.01;
    total_time = 0.1;
    dt = r*dx*dx;
    gridpoints = int((end_x - start_x)/dx - 2);
    timesteps = int(total_time/dt - 1);

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

    for (int i = 0; i < gridpoints; i++) x[i] = dx*(i+1);
    for (int i = 0; i < timesteps; i++) t[i] = dt*i;
    //Initial condition
    for (int i = 0; i < gridpoints; i++) v[0][i] = -x[i];


    if (method == "explicit"){
      Explicit_scheme(v, x, r, gridpoints, timesteps);

      for (int m = 0; m < timesteps; m++){
        for (int j = 0; j < gridpoints; j++){
          v[m][j] += x[j];
        }
      }


      cout << "t[10] = " << t[100] << endl;

      ofile.open(outfilename);
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v[100][i] << endl;
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
<<<<<<< HEAD
          cout << y[i] << endl;
        }


        *y = *v[m];

=======
        }

>>>>>>> 280a0c3419352a54dbfe711b5ed4fa7c58595d2c
        Forward_substitution(a, b, c, y, gridpoints);
        Back_substitution(v[m+1], b, c, y, gridpoints);
      }

      for (int m = 0; m < timesteps; m++){
        for (int j = 0; j < gridpoints; j++){
          v[m][j] += x[j];
        }
      }


      cout << "t[10] = " << t[10] << endl;

      ofile.open(outfilename);
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v[10][i] << endl;
      }

      ofile.close();

    }


    if (method == "CN"){
      double *a, *b, *c;
      double alpha, beta, gamma;

      double *y;

      a = new double[gridpoints];
      b = new double[gridpoints];
      c = new double[gridpoints];
      y = new double[gridpoints];

      alpha = r; beta = 1 - 2*r; gamma = r;


      for (int m = 0; m < timesteps - 1; m++){
        for (int i = 0; i < gridpoints; i++){
            a[i] = -r;
            b[i] = 1.0 + 2*r;
            c[i] = -r;

            if (i == 0){
              y[i] =  beta*v[m][i] + gamma*v[m][i+1];

            } else if (i == gridpoints-1){
              y[i] =  beta*v[m][i] + alpha*v[m][i-1];

            } else {

              y[i] = beta*v[m][i] + gamma*v[m][i+1] + alpha*v[m][i-1];

            }

        }

      Forward_substitution(a, b, c, y, gridpoints);
      Back_substitution(v[m+1], b, c, y, gridpoints);
      }


      for (int m = 0; m < timesteps; m++){
        for (int j = 0; j < gridpoints; j++){
          v[m][j] += x[j];
        }
      }


      cout << "t[10] = " << t[1000] << endl;

      ofile.open(outfilename);
      for (int i = 0; i < gridpoints; i++){
        ofile << x[i] << " " << v[1000][i] << endl;
      }

      ofile.close();

    }


  }

  

  return 0;
}
