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

using namespace std;

ofstream ofile;

double CoulombRepulsionMC_ImportanceSampling(double*);
double CoulombRepulsionMC(double*);


int main(int nargs, char* args[]){
  string outfilename, integration_method;

  outfilename = string(args[1]);
  integration_method = string(args[2]);

  if (integration_method == "1"){
    //Declaration of variables
    int n, local_n, numprocs, my_rank, d;
    double *x;
    double total_integral, local_integral, alpha;
    double total_sigma, local_sigma, local_variance, total_variance, std_mean;
    double time_start, time_end, total_time;
    double a, b;
    double jacobidet, exact, relative_error, func_value;

    //Sets up the uniform distribution for x in [0,1]
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0,1);

    n = atoi(args[3]);
    a = atoi(args[4]);
    b = atoi(args[5]);

    //Hardcode specific parameters
    d = 6;                                                //dimensions    (three of the integrals are trivially solved analytically)
    jacobidet = pow(b-a, d);                              //Jacobideterminant
    x = new double[d];                                    //Integration coordinates x1,...,xd.
    total_integral = 0.;
    local_integral = 0.;
    local_sigma = 0.;
    total_sigma = 0.;
    local_variance = 0.;
    total_variance = 0.;
    std_mean = 0.;
    exact = 5*M_PI*M_PI/(16*16);
    alpha = 4;

    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    local_n = n/numprocs;
    time_start = MPI_Wtime();

    //Ah shit, here we go again... summing up the function values for n_local monte carlo samples locally
    for (int i = 0; i < local_n; i++){
      //Collect a sample of the stochastic varable X = (X1,...,Xd)
      for (int j = 0; j < d; j++){
        //Use the mapping mu(x[j]) = a + (b-a)*x[j], fill
        x[j] = RandomNumberGenerator(gen);                //Random number from uniform distribution
        x[j] = a + (b-a)*x[j];                            //Change of coordinates.
      }
      func_value = CoulombRepulsionMC(x);                 //Compute function value for the sample of X
      local_integral += func_value;                             //computed the contribution to the integral
      local_sigma += func_value*func_value;                     //computes the contribution to the variance
    }

    MPI_Reduce(&local_integral, &total_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sigma, &total_sigma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    total_integral /= ((double) n);
    total_sigma /= ((double) n);
    total_variance = total_sigma - total_integral*total_integral;
    total_integral *= jacobidet;
    std_mean = jacobidet*sqrt(total_variance/( (double) n));
    time_end = MPI_Wtime();
    total_time = time_end - time_start;
    relative_error = abs((total_integral - exact)/exact);

    if (my_rank == 0){
      ofile.open(outfilename);
      ofile << total_integral << " " << std_mean << " " << relative_error << " " << total_time << endl;
      ofile.close();
    }
    MPI_Finalize();
  }

  if (integration_method == "2"){
    //Declaration of variables
    int n, local_n, numprocs, my_rank, d;
    double *x;
    double total_integral, local_integral, alpha;
    double total_sigma, local_sigma, total_variance, local_variance, std_mean;
    double time_start, time_end, total_time, max_radial_distance;
    double jacobidet, exact, relative_error, func_value;

    //Sets up the uniform distribution for x in [0,1]
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0,1);

    n = atoi(args[3]);
    max_radial_distance = atof(args[4]);


    //Hardcode specific parameters
    d = 3;                                                //dimensions    (three of the integrals are trivially solved analytically)
    jacobidet = 8*pow(M_PI, 3)/16;                        //"Jacobideterminant"
    x = new double[d];                                    //Integration coordinates x1,...,xd.
    total_integral = 0.;
    local_integral = 0.;
    local_sigma = 0.;
    total_sigma = 0.;
    local_variance = 0.;
    total_variance = 0.;
    std_mean = 0.;
    exact = 5*M_PI*M_PI/(16*16);
    alpha = 4;



    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    local_n = n/numprocs;
    time_start = MPI_Wtime();

    //Again, it begins... summing up the function values for n_local monte carlo samples locally
    for (int i = 0; i < local_n; i++){
      //Sample from the respective probability distributions of X = (r1,r2,theta2)
      x[0] = RandomNumberGenerator(gen);
      x[0] = -log(1-x[0])/alpha;                      //r1-coordinate
      x[1] = RandomNumberGenerator(gen);
      x[1] = -log(1-x[1])/alpha;                      //r2-coordinate
      x[2] = RandomNumberGenerator(gen);
      x[2] = M_PI*x[2];                               //theta2-coordinate
      func_value = CoulombRepulsionMC_ImportanceSampling(x);
      local_integral += func_value;
      local_sigma += func_value*func_value;
    }


    MPI_Reduce(&local_integral, &total_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sigma, &total_sigma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    total_integral /= ((double) n);                             //Arithmetic mean to obtain the integral
    total_sigma /= ((double) n);                                //Arithmetic mean to obtain the variance
    total_variance = total_sigma - total_integral*total_integral;                 //Computes the variance
    total_integral *= jacobidet;                                //Final contribution to the integral.
    std_mean = jacobidet*sqrt(total_variance/( (double) n));
    time_end = MPI_Wtime();
    total_time = time_end - time_start;
    relative_error = abs((total_integral-exact)/exact);

    //Write to file
    if (my_rank == 0){
      ofile.open(outfilename);
      ofile << total_integral << " " << std_mean << " " << relative_error << " " << total_time << endl;
      ofile.close();
    }

    MPI_Finalize();
  }
  return 0;
}


double CoulombRepulsionMC_ImportanceSampling(double *x){
  double func_value;
  double norm = sqrt( x[0]*x[0] + x[1]*x[1] - 2*x[0]*x[1]*cos(x[2]) );
  if (norm == 1){
    func_value = 0;
  }
  else{
    func_value = x[0]*x[0]*x[1]*x[1]*sin(x[2])/norm;
  }
  return func_value;
}

double CoulombRepulsionMC(double *x){
  double func_value = 0;
  double norm = sqrt((x[0]-x[3])*(x[0]-x[3]) + (x[1]-x[4])*((x[1]-x[4])) + (x[2]-x[5])*(x[2]-x[5]));
  if (norm == 0){
    func_value = 0;
  }
  else{
    func_value = exp(-4*( sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) + sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]) ))/norm;
  }
  return func_value;
}
