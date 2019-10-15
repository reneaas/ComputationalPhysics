#include <iostream>
#include "lib.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <armadillo>
#include <string>
#include "time.h"
#include <random>
#include <mpi.h>

ofstream ofile;

using namespace std;

double CoulombRepulsion(double, double, double, double, double, double);
double CoulombRepulsion_spherical(double, double, double);
double F(double, double);
double radial_probability_density(double,double);
double test_func(double, double, double);
double gauss_legendre_solver(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int n, double (*)(double, double, double, double, double, double));

int main(int nargs, char* args[]){
  //Arguments taken from the command line
  int N = atoi(args[1]);                        //Number of monte carlo samples
  string write_to_file = string(args[2]);       //Write to file or not.
  string outfilename;
  if (write_to_file == "write_to_file"){
    outfilename = string(args[3]);
  }

  int n, local_N, numprocs, my_rank;
  double *a, *b;
  double total_integral, local_integral;
  double time_start, time_end, total_time;
  double max_radial_distance = 10;
  n = 10;
  int d = 6;          //d-dimensional integral.
  double r1, r2, theta2;
  a = new double[d];
  b = new double[d];
  double alpha = 4;

  //r1 endpoints
  a[0] = 0;
  b[0] = max_radial_distance;
  //r2 endpoints
  a[1] = 0;
  b[1] = max_radial_distance;
  //theta1 endpoints
  a[2] = 0;
  b[2] = M_PI;
  //theta2 endpoints
  a[3] = 0;
  b[3] = M_PI;
  //phi1 endpoints
  a[4] = 0;
  b[4] = 2*M_PI;
  //phi2 endpoints
  a[5] = 0;
  b[5] = 2*M_PI;
  double *MC_integrals;
  local_integral = 0.0;
  double relative_error;
  double exact = 5*pow(M_PI,2)/(16*16);

  //Initialize the seed and call the Mersienne algorithm
  random_device rd;
  mt19937_64 gen(rd());
  //Sets up the uniform distribution for x in [0,1]
  uniform_real_distribution<double> RandomNumberGenerator(0,1);
  //Creates the variables


  MPI_Init(&nargs, &args);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  time_start = MPI_Wtime();

  local_N = N/numprocs;
  MC_integrals = new double[local_N];

  for (int i = 0; i < local_N; i++){
    /*
    if (my_rank == 0){
      cout << "Computing for sample = " << numprocs*i << endl;
    }
    */

    for (int j = 0; j < n; j++){
      r1 = RandomNumberGenerator(gen);
      r1 = -log(1-r1)/alpha;
      for (int k = 0; k < n; k++){
        r2 = RandomNumberGenerator(gen);
        r2 = -log(1-r2)/alpha;
        for (int l = 0; l < n; l++){
          theta2 = RandomNumberGenerator(gen);
          theta2 = a[3] + (b[3]-a[3])*theta2;
          MC_integrals[i] += CoulombRepulsion_spherical(r1,r2,theta2)/radial_probability_density(r1,r2);
        }
      }
    }
    local_integral += MC_integrals[i] / (pow( (double) n, (double) 3));
  }

  MPI_Reduce(&local_integral, &total_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total_integral /= (double) N;
  total_integral *=  8*M_PI*M_PI*M_PI;                               //Multiplying by factors due to integration with respect to phi1, phi2 and theta1. Integrand not explicitly dependent on them.
  time_end = MPI_Wtime();
  total_time = time_end - time_start;
  relative_error = abs((total_integral-exact)/exact);
  /*
  if (my_rank == 0){
    cout << "Computed integral = " << total_integral << endl;
    cout << "Analytical value = " << 5*pow(M_PI,2)/(16*16) << endl;
    cout << "time used = " << total_time << endl;
  }
  */
  if (write_to_file == "write_to_file" && my_rank == 0){
    ofile.open(outfilename);
    ofile << N << " " << total_time << " " << total_integral << " " << relative_error << endl;
    ofile.close();
  }

  MPI_Finalize();
  return 0;
}



double CoulombRepulsion(double x1, double y1, double z1, double x2, double y2, double z2){
  /*
  The function which we integrate.
  */
  if (sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) == 0){
    return 0;
  }
  return exp( -4*( sqrt(x1*x1 + y1*y1 + z1*z1) + sqrt(x2*x2 + y2*y2 + z2*z2) ))/sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

double CoulombRepulsion_spherical(double r1, double r2, double theta2){
  /*
  CoulombRepulsion in spherical coordinates. Includes the Jacobi determinant r1^2r2^2 sin(theta2).
  The other coordinates are just multiplied as constants since the integral is not explicitly dependent upon phi1, phi2 or theta1.
  */
  if (sqrt(r1*r1 + r2*r2 - r1*r2*cos(theta2)) == 0){
    return 0;
  }
  return r1*r1*r2*r2*sin(theta2)*exp(-4*(r1+r2)) / ( sqrt(r1*r1 + r2*r2 -2*r1*r2*cos(theta2)) );
}


double radial_probability_density(double r1,double r2){
  /*
  Normalized radial probability distribution.
  */
  return 16*exp(-4*(r1+r2));
}

double test_func(double x, double y, double z){
  return x*x*y*y*z*z*exp(-(x+y+z));
}
