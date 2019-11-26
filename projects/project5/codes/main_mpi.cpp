#include <cmath>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include "time.h"
#include "functions.h"
#include <mpi.h>

using namespace std;
ofstream ofile;


int main(int nargs, char* args[]){

  double **total_v, **local_v, *y;
  int my_rank, numprocs;
  int timesteps, gridpoints, local_h;
  double r, dt, h, total_time;
  double start_x, end_x;
  string outfilename;
  double a, b, local_start,local_end;


  a = 0.; b = 1;

  //Command line arguments
/*
  h = atof(args[2]);
  outfilename = string(args[3]);
*/

  h = 0.01;
  outfilename = "test.txt";


  //Hardcode variables.
  start_x = 0.;
  end_x = 1.;
  r = 0.25;
  total_time = 0.0075;
  dt = r*h*h;
  gridpoints = int((end_x - start_x)/h-2);
  timesteps = int(total_time/dt - 1);

  //Initiate empty solution matrix.
  total_v = new double*[gridpoints];
  for (int j = 0; j < gridpoints; j++){
    total_v[j] = new double[gridpoints];
  }


    for (int i = 0; i < gridpoints; i++){
      for (int j = 0; j < gridpoints; j++){
        total_v[i][j] = 0.;
      }
    }

  y = new double[gridpoints];

  for (int i = 0; i < gridpoints; i++){
    y[i] = h*(i+1);                                             //Position array y-direction
  }


  //Initial condition
  for (int i = 1; i < gridpoints-1; i++){
    for (int j = 1; j < gridpoints-1; j++){
      total_v[i][j] = (a-b)*y[j] - a;
    }
  }

  // MPI initializations
  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  local_v = new double*[gridpoints];
  for (int i = 0; i < gridpoints; i++) local_v[i] = new double[gridpoints];

  local_h = gridpoints/numprocs;

  local_start = 1 + my_rank*local_h;
  local_end = (1 + my_rank)*local_h-1;

  for (int t = 0; t < total_time; t++){
    for (int i = local_start; i < local_end; i++){
      for (int j = 1; j < gridpoints-1; j++){
        local_v[i][j] = (1-4*r)*local_v[i][j] + r*(local_v[i+1][j] + local_v[i-1][j] + local_v[i][j+1] + local_v[i][j-1]);
      }
    }

  //MPI_Gather(&local_v, &total_v, gridpoints*gridpoints, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //MPI_Gather(&local_v, gridpoints*gridpoints, MPI_DOUBLE, total_v, gridpoints*gridpoints, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  /*
  for (int i = 0; i < gridpoints; i++){
    MPI_Reduce(&local_v[i], &total_v[i], gridpoints, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_v, &total_v, gridpoints, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  */
  }
    for (int i = 0; i < gridpoints; i++){
      for (int j = 0; j < gridpoints; j++){
        total_v[i][j] += (b-a)*y[j] + a;
      }
    }

  /*
  if (my_rank == 0){
    ofile.open(outfilename);
    ofile << total_time << " " << h << endl;
    for (int i = 0; i < gridpoints; i++){
      for (int j = 0; j < gridpoints; j++){
      ofile << total_v[i][j] << " ";
      }
      ofile << " " << endl;
    }
    ofile.close();
  }
  */

  MPI_Finalize();
  return 0;
}
