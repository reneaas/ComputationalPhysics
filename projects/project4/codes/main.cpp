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
void initialize_ordered(int, int **, double&, double&);
void initialize_random(int, int **, double&, double&);
void Monte_Carlo_Metropolis(int, int, int **, int, double&, double&, double&, double&, double*, double*, double*, double);

inline int periodic(int coordinate, int dimensions, int step) {
  return (coordinate+dimensions+step) % (dimensions);
}

//Main program
int main(int nargs, char* args[]){

  string outfilename, initialize;
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

  initialize = atos(args[5]);

  //filling spin matrix with s = +1 for all elements
  if (initialize == "ordered"){
  initialize_ordered(n, spin_matrix, E, M);
  }

  //filling spin matrix with arbitrary spin values
  if (initialize == "random"){
    initialize_random(n, spin_matrix, E, M);
  }

  double E_squared, M_squared;                //Change in energy and magnetization

  E_squared = E*E;
  M_squared = M*M;
  J = 1;                                     //Coupling constant

  if (number_of_temperatures == 1){
    double T = atof(args[6]);
    double* expectation_values;
    double* analytical_values;
    double* relative_error;
    expectation_values = new double[7];         //expectation_values = (E, E^2, |M|, |M|^2, M, M^2).
    analytical_values = new double[8];
    relative_error = new double[9];
    double n_spins = (double) n*n;

    //Hardcode initial expectation values to zero.
    for (int i = 0; i < 6; i++){
      expectation_values[i] = 0.0;
    }

    for (int i = 0; i < 6; i++){
      analytical_values[i] = 0.0;
    }

    for (int i = 0; i < 8; i++){
      relative_error[i] = 0.0;
    }

    //Computing the boltzmann distribution for 5 values of dE
    double beta = 1/(T);                //k_B = 1
    for (int i = -8; i < 9; i+=4){
      boltzmann_distribution[i + 8] = exp(-beta*i);
    }

    Monte_Carlo_Metropolis(MC_cycles, n, spin_matrix, J, E, M, E_squared,  M_squared, boltzmann_distribution, expectation_values, analytical_values, beta);

    //Prints exact and computed values to screen for the case T = 1 with a (2 x 2)-lattice.

    cout << "----------------Exact values--------------------------- " << endl;
    double Z = 4*(3 + cosh(8));
    cout << "E = " << -32*sinh(8)/Z << endl;
    cout << " E^2 = " << 256*cosh(8)/Z << endl;
    cout << "|M| = " << 8*(exp(8) + 2)/Z << endl;
    cout << "|M|^2 = " << (32*exp(8) + 4)/Z << endl;
    cout << "M = "<< 0 << endl;
    cout << "M^2 = " << 32*(exp(8) + 1)/Z << endl;
    cout << "Heat capacity = " << (256*cosh(8)/Z - (-32*sinh(8)/Z * -32*sinh(8)/Z))/(T*T) << endl;
    cout << "Magnetic susceptibility = " << 32*(exp(8) + 1)/Z/T << endl;
    cout << "-------------------Computed Values -------------------------" << endl;
    cout << "E = " << expectation_values[0] << endl;              // E =
    cout << "E^2 = " << expectation_values[1] << endl;
    cout << "|M| = " << expectation_values[2] << endl;
    cout << "|M|^2 = " << expectation_values[3] << endl;
    cout << "M = " << expectation_values[4] << endl;
    cout << "M^2 = " << expectation_values[5] << endl;
    cout << "Heat capacity = " << (expectation_values[1] - (expectation_values[0]*expectation_values[0]))/(T*T) << endl;
    cout << "Magnetic susceptibility = " << expectation_values[5]/T << endl;



  }


  if (number_of_temperatures > 1){
    T_initial = atof(args[6]);                      //Initial temperature
    T_final = atof(args[7]);                        //Final temperature
    step_size = atof(args[8]);                      //temperature step size.

    double** expectation_values;                    //matrix to store computed expectation values.
    double** analytical_values;                     //matrix to store analytical expectation values.
    double** relative_error;                        //matrix to store relative error.
    int** initial_spin_matrix;                      //Stores the initial spin matrix.
    double E_initial, M_initial;                    //Stores initial energy and magnetization of system.
    double* magnetic_susceptibility;                //Stores the computed magnetic susceptibilities for each temperature
    double* magnetic_susceptibility_analytical;     //Stores the analytical magnetic susceptibilities for each temperature
    double* heat_capacity;                          //Stores the computed heat capacity for each temperature.
    double* heat_capacity_analytical;               //Stores the analytical heat capacity for each temperature.
    int n_spins = n*n;                              //Total number of spins.

    magnetic_susceptibility = new double[number_of_temperatures + 1];
    heat_capacity = new double[number_of_temperatures + 1];

    magnetic_susceptibility_analytical = new double[number_of_temperatures + 1];
    heat_capacity_analytical = new double[number_of_temperatures + 1];


    //Store the initial values of the system.
    E_initial = E;
    M_initial = M;
    initial_spin_matrix = spin_matrix;

    expectation_values = new double*[number_of_temperatures + 1];       //a vector for computed expectation values for each temperature.
    for (int i = 0; i <= number_of_temperatures; i++){
      expectation_values[i] = new double[5];
    }

    analytical_values = new double*[number_of_temperatures + 1];        //a vector for analytical expectation values for each temperature.
    for (int i = 0; i <= number_of_temperatures; i++){
      analytical_values[i] = new double[5];
    }

    relative_error = new double*[number_of_temperatures + 1];        //a vector for analytical expectation values for each temperature.
    for (int i = 0; i <= number_of_temperatures; i++){
      relative_error[i] = new double[7];
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
      Monte_Carlo_Metropolis(MC_cycles, n, spin_matrix, J, E, M, E_squared,  M_squared, boltzmann_distribution, expectation_values[i], analytical_values[i],beta);

      //Computing magnetic susceptibility and heat capacity for each temperature
      heat_capacity_analytical[i] = (analytical_values[i][1]-analytical_values[i][0]*analytical_values[i][0])*beta*beta;     //Stores the analytical expectation value for heat capacity for a given temperature
      magnetic_susceptibility_analytical[i] = (analytical_values[i][5])*beta;                                                //Stores the analytical expectation value for susceptibility for a given temperature


      heat_capacity[i] = (expectation_values[i][1]-expectation_values[i][0]*expectation_values[i][0])*beta*beta;             //Stores the computed expectation value for heat capacity for a given temperature
      magnetic_susceptibility[i] = (expectation_values[i][5]-(expectation_values[i][4]*expectation_values[i][4]))*beta;        //Stores the computed expectation value for susceptibility for a given temperature

      relative_error[i][0] = abs((analytical_values[i][0]-expectation_values[i][0])/analytical_values[i][0]);    //Stores relative error in E
      relative_error[i][1] = abs((analytical_values[i][1]-expectation_values[i][1])/analytical_values[i][1]);    //Stores relative error in E_squared
      relative_error[i][2] = abs((analytical_values[i][2]-expectation_values[i][2])/analytical_values[i][2]);    //Stores relative error in M_abs
      relative_error[i][3] = abs((analytical_values[i][3]-expectation_values[i][3])/analytical_values[i][3]);    //Stores relative error in M_abs_squared
      relative_error[i][4] = abs((analytical_values[i][4]-expectation_values[i][4])/analytical_values[i][4]);    //Stores relative error in M
      relative_error[i][5] = abs((analytical_values[i][5]-expectation_values[i][5])/analytical_values[i][5]);    //Stores relative error in M_squared
      relative_error[i][6] = abs((heat_capacity_analytical[i]-heat_capacity[i])/heat_capacity_analytical[i]);    //Stores relative error in Heat Capacity
      relative_error[i][7] = abs((magnetic_susceptibility_analytical[i]-magnetic_susceptibility[i])/magnetic_susceptibility_analytical[i]);  //Stores relative error in Heat Capacity


      //Prints exact and computed values to screen
      cout << "----------------T = " << temperatures[i] <<"---------------" << endl;
      cout << "----------------Magnetic susceptibility------------- " << endl;
      cout << "Computed = " << magnetic_susceptibility[i] << endl;
      cout << "Analytical = " << magnetic_susceptibility_analytical[i] << endl;
      cout << "-------------------Heat capacity ------------------" << endl;
      cout << "Computed = " << heat_capacity[i] << endl;
      cout << "Analytical = " << heat_capacity_analytical[i] << endl;


    }

    //Printing out just to test.
    for (int i = 0; i <= number_of_temperatures; i++){
      cout << "T = " << temperatures[i] << " ," << "E = " << expectation_values[i][0]/((double) n_spins) << endl;
    }


    //Write the computed expectation values to file here...

  }

  return 0;
}

void initialize_ordered(int dimensions, int **spin_matrix, double& E, double& M){
  // setup spin matrix and intial magnetization
  for(int i =0; i <dimensions; i++) {
    for (int j= 0; j < dimensions; j++){
      spin_matrix[i][j] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[i][j];
    }
  }
  // setup initial energy
  for(int i = 0; i < dimensions; i++) {
    for (int j= 0; j < dimensions; j++){
      E -=  (double) spin_matrix[i][j] * (spin_matrix[periodic(i,dimensions,-1)][j] + spin_matrix[i][periodic(j,dimensions,-1)]);
      }
    }
  }


void initialize_random(int dimensions, int **spin_matrix, double& E, double& M){
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0,1);

    double x;

    // setup spin matrix and intial magnetization
    for(int i = 0; i <dimensions; i++) {
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
    // setup initial energy
    for(int i =0; i < dimensions; i++) {
      for (int j= 0; j < dimensions; j++){
        E -=  (double) spin_matrix[i][j] * (spin_matrix[periodic(i,dimensions,-1)][j] + spin_matrix[i][periodic(j,dimensions,-1)]);
        }
      }
  }


void Monte_Carlo_Metropolis(int MC, int n, int **spin_matrix, int J, double& E, double& M, double& E_squared, double& M_squared,
                            double* boltzmann_distribution, double* expectation_values, double* analytical_values, double beta){

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
  expectation_values[5] = M_squared;

  double Z_a,E_a,M_a,E_squared_a,M_squared_a, Mabs_a, Mabs_squared_a;

  Z_a =  4*(3 + cosh(8*beta));
  E_a = -32*sinh(8*beta)/Z_a;
  E_squared_a = 256*cosh(8*beta)/Z_a;
  M_a = 0;
  M_squared_a = 32*(exp(8*beta) + 1)*beta/Z_a ;
  Mabs_a =  8*(exp(8*beta) + 2)/Z_a;
  Mabs_squared_a = (32*exp(8*beta) + 4)/Z_a;

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

}
