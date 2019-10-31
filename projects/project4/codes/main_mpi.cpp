#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <string>
#include "time.h"
#include <random>
#include <fstream>
#include <mpi.h>
using namespace  std;

inline int periodic(int coordinate, int dimensions, int step) {
  return (coordinate+dimensions+step) % (dimensions);

int main(int nargs, char* args[]){
  string outfilename;
  double boltzmann_distribution[17];
  int n, MC_samples, , J;
  int **spin_matrix, my_rank, numprocs;
  double E_initial, M_initial;                    //Stores initial energy and magnetization of system.
  int n_spins;                                     //Total number of spins.
  double beta;
  double*


  //Read from command line
  outfilename = string(args[1]);         //Name of the file to write the results to
  n = atoi(args[2]);                    //Dimension of spin matrix
  MC_samples = atoi(args[3]);            //Number of Monte Carlo samples


  //initialize matrix
  spin_matrix = new int*[n];
  for (int i = 0; i < n; i++){
    spin_matrix[i] = new int[n];
  }

  E_initial = 0;
  M_initial = 0;
  n_spins = n*n;
  J = 1;                                     //Coupling constant
  initialize(n, spin_matrix, E_initial, M_initial, "random");


// MPI initializations
MPI_Init (&nargs, &args);
MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  return 0;

}



void initialize(int dimensions, int **spin_matrix, double& E, double& M, string initialize){
  if (initialize == "ordered"){
    // setup spin matrix and intial magnetization
    for(int i =0; i < dimensions; i++) {
      for (int j= 0; j < dimensions; j++){
        spin_matrix[i][j] = 1; // spin orientation for the ground state
        M +=  (double) spin_matrix[i][j];
      }
    }
  }

  if (initialize == "random"){
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0,1);
    double x;

    // setup spin matrix and intial magnetization
    for(int i = 0; i < dimensions; i++) {
      for (int j= 0; j < dimensions; j++){
        x = RandomNumberGenerator(gen);
        if (x <= 0.5){
          spin_matrix[i][j] = 1;        //spin orientation for the ground state
        }
        else{
          spin_matrix[i][j] = -1;      //spin orientation for the ground state
        }
        M +=  (double) spin_matrix[i][j];
      }
    }

  }
  // setup initial energy
  for(int i = 0; i < dimensions; i++) {
    for (int j= 0; j < dimensions; j++){
      E -=  (double) spin_matrix[i][j] * (spin_matrix[periodic(i,dimensions,-1)][j] + spin_matrix[i][periodic(j,dimensions,-1)]);
      }
    }
  }
