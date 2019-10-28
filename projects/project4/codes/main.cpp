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

ofstream ofile;               //Global variable for writing results to file.


//Declaration of functions.
void initialize(int, int **, double&, double&);
void Monte_Carlo_Metropolis(int, int, int **, int, double&, double&, double&, double&, double*, double*);

inline int periodic(int coordinate, int dimensions, int step) {
  return (coordinate+dimensions+step) % (dimensions);
}

//Main program
int main(int nargs, char* args[]){

<<<<<<< HEAD
  string outfilename, number_of_temperatures;
  double E, M,boltzmann_distribution[9];
=======
  string outfilename;
  double E, M,boltzmann_distribution[17];
>>>>>>> f120a5b26d26e3bdb94c414b49193a8286232487
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

<<<<<<< HEAD
  //Computing the boltzmann distribution for 5 values of dE
  double beta = 1/(T);                //k_B = 1
  for (int i = 0; i < 9; i+=4){
    boltzmann_distribution[i] = exp(-beta*i);
  }

=======
>>>>>>> f120a5b26d26e3bdb94c414b49193a8286232487
  E = 0;
  M = 0;

  //filling spin matrix with arbitrary spin values
  initialize(n, spin_matrix, E, M);


<<<<<<< HEAD
    int x_flip, y_flip, J;
    double dE, dM, dE_squared, dM_squared;                //Change in energy and magnetization

    dE_squared = 0;
    dM_squared = 0;

    J = 1;
=======

  double E_squared, M_squared;                //Change in energy and magnetization
>>>>>>> f120a5b26d26e3bdb94c414b49193a8286232487

  E_squared = E*E;
  M_squared = M*M;
  J = 1;                                     //Coupling constant

  if (number_of_temperatures == 1){
    double T = atof(args[5]);
    double* expectation_values;
    expectation_values = new double[6];         //expectation_values = (E, E^2, |M|, |M|^2, M, M^2).
    double n_spins = (double) n*n;

<<<<<<< HEAD
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
=======
    //Hardcode initial expectation values to zero.
    for (int i = 0; i < 6; i++){
      expectation_values[i] = 0.0;
    }
    //Computing the boltzmann distribution for 5 values of dE
    double beta = 1/(T);                //k_B = 1
    for (int i = -8; i < 9; i+=4){
      boltzmann_distribution[i + 8] = exp(-beta*i);
    }

    Monte_Carlo_Metropolis(MC_cycles, n, spin_matrix, J, E, M, E_squared,  M_squared, boltzmann_distribution, expectation_values);

    //Prints exact and computed values to screen for the case T = 1 with a (2 x 2)-lattice.
    cout << "----------------Exact values--------------------------- " << endl;
    double Z = 4*(3 + cosh(8));
    cout << "E = " << -32*sinh(8)/Z << endl;
    cout << " E^2 = " << 256*cosh(8)/Z << endl;
    cout << "|M| = " << 8*(exp(8) + 2)/Z << endl;
    cout << "|M|^2 = " << (32*exp(8) + 4)/Z << endl;
    cout << "M = "<< 0 << endl;
    cout << "M^2 = " << 32*(exp(8) + 1)/Z << endl;
    cout << "-------------------Computed Values -------------------------" << endl;
    cout << "E = " << expectation_values[0] << endl;              // E =
    cout << "E^2 = " << expectation_values[1] << endl;
    cout << "|M| = " << expectation_values[2] << endl;
    cout << "|M|^2 = " << expectation_values[3] << endl;
    cout << "M = " << expectation_values[4] << endl;
    cout << "M^2 = " << expectation_values[5] << endl;

>>>>>>> f120a5b26d26e3bdb94c414b49193a8286232487

  }


  if (number_of_temperatures > 1){
    T_initial = atof(args[5]);          //Initial temperature
    T_final = atof(args[6]);            //Final temperature
    step_size = atof(args[7]);          //temperature step size.
    double** expectation_values;        //matrix to store expectation values.
    int** initial_spin_matrix;          //Stores the initial spin matrix.
    double E_initial, M_initial;        //Stores initial energy and magnetization of system.
    double* magnetic_susceptibility;    //Stores the computed magnetic susceptibilities for each temperature
    double* heat_capacity;              //Stores heat capacity for each temperature.
    int n_spins = n*n;                  //Total number of spins.

    magnetic_susceptibility = new double[number_of_temperatures + 1];
    heat_capacity = new double[number_of_temperatures + 1];


    //Store the initial values of the system.
    E_initial = E;
    M_initial = M;
    initial_spin_matrix = spin_matrix;

    expectation_values = new double*[number_of_temperatures + 1];       //a vector for computed expectation values for each temperature.
    for (int i = 0; i <= number_of_temperatures; i++){
      expectation_values[i] = new double[5];
    }

    analytical_values = new double*[number_of_temperatures + 1];
    for (int i = 0; i <= number_of_temperatures; i++){
      analytical_values[i] = new double[5];
    }



    //Store temperatures to loop over.
    double* temperatures;
    temperatures = new double[number_of_temperatures+1];
    for (int i = 0; i <= number_of_temperatures; i++){
      temperatures[i] = T_initial + (double) step_size*i;
    }

    for (int i = 0; i <= number_of_temperatures; i++){
      //Makes sure that we use the same initial state for each simulation.
      spin_matrix = initial_spin_matrix;
      E = E_initial;
      M = M_initial;

      //Computing the boltzmann distribution for 5 values of dE
      double beta = 1/((double) temperatures[i]);                //k_B = 1
      for (int j = -8; j < 9; j+=4){
        boltzmann_distribution[j + 8] = exp(-beta*j);
      }

      //Metro time.
      Monte_Carlo_Metropolis(MC_cycles, n, spin_matrix, J, E, M, E_squared,  M_squared, boltzmann_distribution, expectation_values[i]);

      //Compute magnetic susceptibility and heat capacity for each temperature here...

    }
<<<<<<< HEAD

    E += dE;                //Computing change in energy
    M += dM;                //Computing change in magnetization
    dE_squared += dE*dE;
    dM_squared += dM*dM;
  }
=======
>>>>>>> f120a5b26d26e3bdb94c414b49193a8286232487

    //Printing out just to test.
    for (int i = 0; i <= number_of_temperatures; i++){
      cout << "T = " << temperatures[i] << " ," << "E = " << expectation_values[i][0]/((double) n_spins) << endl;
    }


    //Write the computed expectation values to file here...

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


<<<<<<< HEAD
void Monte_Carlo_Metropolis(){


=======
void Monte_Carlo_Metropolis(int MC, int n, int **spin_matrix, int J, double& E, double& M, double& E_squared, double& M_squared,
                            double* boltzmann_distribution, double* expectation_values){

  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntegerGenerator(0,n-1);        //Sets up the uniform distribution for x in [0,n-1]
  uniform_real_distribution<double> RandomNumberGenerator(0,1);       //Sets up the uniform distribution for x in [0,1]

  int x_flip, y_flip, dE, dM, n_spins;
  double E_sum, M_sum, Mabs_sum, Mabs_sum_squared;

  n_spins = n*n;

  E_sum = 0;
  M_sum = 0;
  Mabs_sum = 0.;
  Mabs_sum_squared = 0.;

  //Running over Monte Carlo samples
  for (int k = 0; k < MC; k++){


    for (int j = 0; j < n_spins; j++){
      x_flip = RandomIntegerGenerator(gen);
      y_flip = RandomIntegerGenerator(gen);      //Randomized indices to matrix element that will be flipped


      dE =  2*J*spin_matrix[x_flip][y_flip] * (spin_matrix[periodic(x_flip,n,1)][y_flip]
            + spin_matrix[periodic(x_flip, n,-1)][y_flip]
            + spin_matrix[x_flip][periodic(y_flip, n,1)]
            + spin_matrix[x_flip][periodic(y_flip, n,-1)]);

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
        dE = 0;                             //Rejecting the flip
        dM = 0;
      }

      E += dE;
      M += dM;

    }

    E_sum += (double) E;
    M_sum += (double) M;
    Mabs_sum += (double) abs(M);

    E_squared += (double) E*E;
    M_squared += (double) M*M;
    Mabs_sum_squared += (double) abs(M)*abs(M);
  }

  E_sum /= (double) (MC);
  M_sum /= (double) MC;
  E_squared /= (double) MC;
  M_squared /= (double) MC;
  Mabs_sum /= (double) MC;
  Mabs_sum_squared /= (double) MC;

  //Store the computed expectation values
  expectation_values[0] = E_sum;
  expectation_values[1] = E_squared;
  expectation_values[2] = Mabs_sum;
  expectation_values[3] = Mabs_sum_squared;
  expectation_values[4] = M_sum;
  expectation_values[5] = Mabs_sum_squared;

  double Z_a,E_a,M_a,E_squared_a,M_squared_a, Mabs_a, Mabs_squared_a;

  Z_a =  4*(3 + cosh(8/T));
  E_a = -32*sinh(8/T)/Z;
  E_squared_a = 256*cosh(8*Tk_inv)/Z
  M_a = 0;
  M_squared_a = 32*(exp(8*Tk_inv) + 1)*Tk_inv/Z ;
  Mabs_a =  8*(exp(8*Tk_inv) + 2)/Z;
  Mabs_squared_a = (32*exp(8*Tk_inv) + 4)/Z);

  analytical_values[0] = E_a;
  analytical_values[1] = E_squared_a;
  analytical_values[2] = Mabs_a;
  analytical_values[3] = Mabs_squared_a;
  analytical_values[4] = M_a;
  analytical_values[5] = M_squared_a;

  /*
  //Store the computed expectation values per spin.
  expectation_values[0] = E_sum/((double) n_spins);
  expectation_values[1] = E_squared/((double) n_spins);
  expectation_values[2] = Mabs_sum/((double) n_spins);
  expectation_values[3] = Mabs_sum_squared/((double) n_spins);
  expectation_values[4] = M_sum/((double) n_spins);
  expectation_values[5] = Mabs_sum_squared/((double) n_spins);
  */
>>>>>>> f120a5b26d26e3bdb94c414b49193a8286232487

}
