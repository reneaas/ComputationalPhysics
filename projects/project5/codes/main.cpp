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
    r = 1;
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
          a[i] = -1.0;
          b[i] = 1.0 + 2*r;
          c[i] = -1.0;
          //cout << y[i] << endl;
          y[i] = v[m][i];
          cout << y[i] << endl;
        }


        *y = *v[m];

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

  }
  return 0;
}
