#include <iostream>
#include "lib.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <string>
#include "time.h"
#include <random>
#include <fstream>
#include <omp.h>
using namespace  std;


void initialize(int, int **, double&, double&);


ofstream ofile;               //Global variable for writing results to file.

inline int periodic(int coordinate, int dimensions, int step) {
  return (coordinate+dimensions+step) % (dimensions);
}

int main(int nargs, char* args[]){

  string outfilename, number_of_temperatures;
  double E, M,boltzmann_distribution[9];
  double T_initial, T_final, step_size;   //Variables to define temperature vector
  int n, MC_cycles;

  //Read from command line
  number_of_temperatures = atoi(args[1]);

  outfilename = string(args[2]);         //Name of the file to write the results to
  n = atoi(args[3]);                    //Dimension of spin matrix
  MC_cycles = atoi(args[4]);            //Number of Monte Carlo cycles

  if (number_of_temperatures != 1){
    T_initial = atof(args[5]);         //Inital temperature
    T_final = atof(args[6]);           //Final temperature
    step_size = atof(args[7]);
  }

  //initialize matrix
  int **spin_matrix;
  spin_matrix = new int*[n];
  for (int i = 0; i < n; i++){
    spin_matrix[i] = new int[n];
  }

  //Computing the boltzmann distribution for 5 values of dE
  double beta = 1/(T);                //k_B = 1
  for (int i = 0; i < 9; i+=4){
    boltzmann_distribution[i] = exp(-beta*i);
  }

  E = 0;
  M = 0;

  //filling spin matrix with arbitrary spin values
  initialize(n, spin_matrix, E, M);

  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntegerGenerator(0,n-1);        //Sets up the uniform distribution for x in [0,n-1]
  uniform_real_distribution<double> RandomNumberGenerator(0,1);       //Sets up the uniform distribution for x in [0,1]

    int x_flip, y_flip, J;
    double dE, dM, dE_squared, dM_squared;                //Change in energy and magnetization

    dE_squared = 0;
    dM_squared = 0;

    J = 1;

  //Running over Monte Carlo samples
  for (int k = 0; k < MC; k++){

    x_flip = RandomIntegerGenerator(gen);
    y_flip = RandomIntegerGenerator(gen);      //Randomized indices to matrix element that will be flipped

    spin_matrix[x_flip][y_flip] *= (-1);

    dM = M + 2*spin_matrix[x_flip][y_flip]*spin_matrix[x_flip][y_flip];

    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        dE = (double) 2*J*spin_matrix[x_flip][y_flip] * (spin_matrix[periodic(x_flip,dimensions,1)][y_flip] + spin_matrix[periodic(x_flip, dimensions,-1)][y_flip] + spin_matrix[x_flip][periodic(y_flip, dimensions,1)] + spin_matrix[x_flip][periodic(y_flip, dimensions,-1)]);
      }
    }

    //Metropolis algorithm
    if (dE >= 0){
      r = RandomNumberGenerator(gen);
      P = boltzmann_distribution[dE];          //Probability from boltzmann distribution
      if (r > P){
        spin_matrix[x_flip][y_flip] *= (-1);   //Rejecting the flip

      }
    }

    E += dE;                //Computing change in energy
    M += dM;                //Computing change in magnetization
    dE_squared += dE*dE;
    dM_squared += dM*dM;
  }


  return 0;
}

void initialize(int dimensions, int **spin_matrix, double& E, double& M){
  // setup spin matrix and intial magnetization
  for(int i =0; i <dimensions; i++) {
    for (int j= 0; j < dimensions; j++){
      spin_matrix[i][j] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[i][j];
    }
  }
  // setup initial energy
  for(int i =0; i < dimensions; i++) {
    for (int j= 0; j < dimensions; j++){
      E -=  (double) spin_matrix[i][j] * (spin_matrix[periodic(i,dimensions,-1)][j] + spin_matrix[i][periodic(j,dimensions,-1)]);
      }
    }
  }


void Monte_Carlo_Metropolis(){



}
