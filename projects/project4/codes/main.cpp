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

void Monte_Carlo_Metropolis(int, int, int **, int, double&, double&, double&, double&, double*);

ofstream ofile;               //Global variable for writing results to file.

inline int periodic(int coordinate, int dimensions, int step) {
  return (coordinate+dimensions+step) % (dimensions);
}

int main(int nargs, char* args[]){

  string outfilename;
  double E, M,boltzmann_distribution[17];
  double T_initial, T_final, step_size;   //Variables to define temperature vector
  int n, MC_cycles, J, number_of_temperatures;

  //Read from command line
  number_of_temperatures = atoi(args[1]);

  outfilename = string(args[2]);         //Name of the file to write the results to
  n = atoi(args[3]);                    //Dimension of spin matrix
  MC_cycles = atoi(args[4]);            //Number of Monte Carlo cycles


  //initialize matrix
  int **spin_matrix;
  spin_matrix = new int*[n];
  for (int i = 0; i < n; i++){
    spin_matrix[i] = new int[n];
  }

  E = 0;
  M = 0;

  //filling spin matrix with arbitrary spin values
  initialize(n, spin_matrix, E, M);



  double E_squared, M_squared;                //Change in energy and magnetization

  E_squared = E*E;
  M_squared = M*M;


  J = 1;                                     //Coupling constant

  if (number_of_temperatures == 1){
    double T = atof(args[5]);

    //Computing the boltzmann distribution for 5 values of dE
    double beta = 1/(T);                //k_B = 1
    for (int i = -8; i < 9; i+=4){
      boltzmann_distribution[i + 8] = exp(-beta*i);
    }

    Monte_Carlo_Metropolis(MC_cycles, n, spin_matrix, J, E, M, E_squared,  M_squared, boltzmann_distribution);


  }


  if (number_of_temperatures > 1){
    T_initial = atof(args[5]);         //Inital temperature
    T_final = atof(args[6]);           //Final temperature
    step_size = atof(args[7]);

    //Store temperatures to loop over.
    double* temperatures;
    temperatures = new double[number_of_temperatures+1];
    for (int i = 0; i <= number_of_temperatures; i++){
      temperatures[i] = T_initial + (double) step_size*i;
    }

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


void Monte_Carlo_Metropolis(int MC, int n, int **spin_matrix, int J, double& E, double& M, double& E_squared, double& M_squared, double* boltzmann_distribution){

  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntegerGenerator(0,n-1);        //Sets up the uniform distribution for x in [0,n-1]
  uniform_real_distribution<double> RandomNumberGenerator(0,1);       //Sets up the uniform distribution for x in [0,1]

  int x_flip, y_flip, dE, dM, nn;
  double E_sum, M_sum;

  nn = n*n;

  E_sum = 0;
  M_sum = 0;
  //Running over Monte Carlo samples
  for (int k = 0; k < MC; k++){


    for (int j = 0; j < nn; j++){
      x_flip = RandomIntegerGenerator(gen);
      y_flip = RandomIntegerGenerator(gen);      //Randomized indices to matrix element that will be flipped


      dE =  2*J*spin_matrix[x_flip][y_flip] * (spin_matrix[periodic(x_flip,n,1)][y_flip] + spin_matrix[periodic(x_flip, n,-1)][y_flip] + spin_matrix[x_flip][periodic(y_flip, n,1)] + spin_matrix[x_flip][periodic(y_flip, n,-1)]);

      //Metropolis algorithm
      if(dE < 0){
        spin_matrix[x_flip][y_flip] *= (-1);   //Accepting the flip
        dM = 2*spin_matrix[x_flip][y_flip];
      }
      else if(RandomNumberGenerator(gen) < boltzmann_distribution[dE + 8]){
        spin_matrix[x_flip][y_flip] *= (-1);     //Accepting flip
        dM = 2*spin_matrix[x_flip][y_flip];
      }
      else{
        dE = 0;
        dM = 0;
      }

      E += dE;
      M += dM;

    }

    E_sum += (double) E;
    M_sum += (double) M;

    E_squared += (double) E*E;
    M_squared += (double) M*M;
  }

  E_sum /= (double) (MC);
  M_sum /= (double) MC;
  E_squared /= MC;
  M_squared /= MC;

  cout<<"Magnetization:" << M_sum <<endl;
  cout << "M_squared = " << M_squared << endl;
  cout << "E = " << E_sum << endl;
  cout << "E^2 = " << E_squared << endl;

}
