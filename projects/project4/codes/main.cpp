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

ofstream ofile, ofile2;               //Global variable for writing results to file.


//Declaration of functions.
void initialize(int, int **, double&, double&, string);
void Monte_Carlo_Metropolis_time(int, int, int **, int, double&, double&, double*, double*, double*, double*, double);
void analytical_values_2x2Lattice(double*, double);
void Monte_Carlo_Metropolis_2x2(int, int, int **, int, double&, double&, double*, double*, double*, double, double**);



inline int periodic(int coordinate, int dimensions, int step) {
  return (coordinate+dimensions+step) % (dimensions);
}

//Main program
int main(int nargs, char* args[]){

  string outfilename, initialize;
  double E, M,boltzmann_distribution[17];
  double T_initial, T_final, step_size;   //Variables to define temperature vector
  int n, MC_samples, J, number_of_temperatures;
  double E_squared, M_squared;
  int **spin_matrix;
  int** initial_spin_matrix;                      //Stores the initial spin matrix.
  double E_initial, M_initial;                    //Stores initial energy and magnetization of system.
  int n_spins;                                     //Total number of spins.
  double beta;

  //Read from command line
  number_of_temperatures = atoi(args[1]);
  outfilename = string(args[2]);         //Name of the file to write the results to
  n = atoi(args[3]);                    //Dimension of spin matrix
  MC_samples = atoi(args[4]);            //Number of Monte Carlo samples
  initializing = string(args[5]);


  //initialize matrix
  initial_spin_matrix = new int*[n];
  for (int i = 0; i < n; i++){
    initial_spin_matrix[i] = new int[n];
  }

  E_initial = 0;
  M_initial = 0;
  n_spins = n*n;
  J = 1;                                     //Coupling constant
  initialize(n, initial_spin_matrix, E_initial, M_initial, initializing);


  if (number_of_temperatures == 1){
    double magnetic_susceptibility;                //Stores the computed magnetic susceptibilities for each temperature
    double heat_capacity;                          //Stores the computed heat capacity for each temperature.
    double *energy, *magnetization, *time, *acceptance;
    double T = atof(args[6]);
    int n_times = MC_samples/n_spins;


    energy = new double[n_times];
    magnetization = new double[n_times];
    time = new double[n_times];
    acceptance = new double[n_times];

    //Compute Boltzmann factors.
    beta = 1/(T);                //k_B = 1
    for (int i = -8; i < 9; i += 4){
      boltzmann_distribution[i + 8] = exp(-beta*i);
    }


    Monte_Carlo_Metropolis_time(MC_samples, n, spin_matrix, J, E, M, boltzmann_distribution, energy, magnetization, time, acceptance, beta);

    ofile.open(outfilename);
    for (int i = 0; i < n_times; i++){
      ofile << time[i] << " " << setprecision(8) << energy[i] << " " << setprecision(8) << " " << magnetization[i] << " " << acceptance[i] << endl;
    }


  }


  if (number_of_temperatures == 1 && number_of_MC_runs == 1){
    double T = atof(args[7]);
    double* expectation_values;
    expectation_values = new double[6];         //expectation_values = (E, E^2, |M|, |M|^2, M, M^2).


    //Hardcode initial expectation values to zero.
    for (int i = 0; i < 6; i++){
      expectation_values[i] = 0.0;
    }

    //Computing the boltzmann distribution for 5 values of dE
    double beta = 1/(T);                //k_B = 1
    for (int i = -8; i < 9; i+=4){
      boltzmann_distribution[i + 8] = exp(-beta*i);
    }

    Monte_Carlo_Metropolis(MC_samples, n, spin_matrix, J, E_initial, M_initial, boltzmann_distribution, expectation_values, beta);


    cout << "-------------------Computed Values-------------------------" << endl;
    cout << "E = " << expectation_values[0] << endl;              // E =
    cout << "E^2 = " << expectation_values[1] << endl;
    cout << "|M| = " << expectation_values[2] << endl;
    cout << "|M|^2 = " << expectation_values[3] << endl;
    cout << "M = " << expectation_values[4] << endl;
    cout << "M^2 = " << expectation_values[5] << endl;



  }


  if (number_of_temperatures > 1 && number_of_MC_runs == 1){
    T_initial = atof(args[7]);                      //Initial temperature
    step_size = atof(args[8]);                      //temperature step size.

    double** expectation_values;                    //matrix to store computed expectation values.
    double* magnetic_susceptibility;                //Stores the computed magnetic susceptibilities for each temperature
    double* heat_capacity;                          //Stores the computed heat capacity for each temperature.

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

    //Store temperatures to loop over.
    double* temperatures;
    temperatures = new double[number_of_temperatures+1];
    for (int i = 0; i <= number_of_temperatures; i++){
      temperatures[i] = T_initial + (double) step_size*i;
    }


    for (int i = 0; i <= number_of_temperatures; i++){
      //Makes sure that we use the same initial state for each simulation.
      spin_matrix =  initial_spin_matrix;
      E = E_initial;
      M = M_initial;

      //Computing the boltzmann distribution for 5 values of dE
      double beta = 1/((double) temperatures[i]);                //k_B = 1
      for (int j = -8; j < 9; j+=4){
        boltzmann_distribution[j + 8] = exp(-beta*j);
      }

      //Metro time.
      Monte_Carlo_Metropolis(MC_samples, n, spin_matrix, J, E, M, boltzmann_distribution, expectation_values[i],beta);

      //Computing magnetic susceptibility and heat capacity for each temperature

      heat_capacity[i] = (expectation_values[i][1]-expectation_values[i][0]*expectation_values[i][0])*beta*beta;             //Stores the computed expectation value for heat capacity for a given temperature
      magnetic_susceptibility[i] = (expectation_values[i][5]-(expectation_values[i][4]*expectation_values[i][4]))*beta;        //Stores the computed expectation value for susceptibility for a given temperature


      //Prints exact and computed values to screen
      cout << "----------------T = " << temperatures[i] <<"---------------" << endl;
      cout << "----------------Magnetic susceptibility------------- " << endl;
      cout << "Computed = " << magnetic_susceptibility[i] << endl;
      cout << "-------------------Heat capacity ------------------" << endl;
      cout << "Computed = " << heat_capacity[i] << endl;

    }

    //Printing out just to test.
    for (int i = 0; i <= number_of_temperatures; i++){
      cout << "T = " << temperatures[i] << " ," << "E = " << expectation_values[i][0]/((double) n_spins) << endl;
    }


    //Write the computed expectation values to file here...

  }


  if (n == 2 && number_of_temperatures == 1){
    double *time, *acceptance;
    double T = atof(args[6]);
    double** expectation_values;
    double*analytical_values;
    double** relative_error;
    int n_times;
    string outfilename = "Expectation_values_n_2.txt"
    string outfilename2 = "Relative_error_n_2.txt";

    energy = new double[n_times];
    magnetization = new double[n_times];
    time = new double[n_times];
    acceptance = new double[n_times];

    n_times = MC_samples/n_spins;
    analytical_values = new double[8];
    expectation_values = new double*[8];
    relative_error = new double*[7];

    for (int i = 0; i < 8; i++){
      expectation_values[i] = new double[n_times];
    }

    for (int i=0; i < 7; i++){
      relative_error[i] = new double[n_times];
    }


    //Computing the boltzmann distribution for 5 values of dE
    double beta = 1/(T);                //k_B = 1
    for (int i = -8; i < 9; i+=4){
      boltzmann_distribution[i + 8] = exp(-beta*i);
    }


    Monte_Carlo_Metropolis_2x2( MC, n, spin_matrix, J,  E, M, boltzmann_distribution, time, acceptance, beta, expectation_values[i]);

    analytical_values_2x2Lattice(analytical_values, T);

    ofile.open(outfilename);
    ofile2.open(outfilename2);
    for (int i = 0; i <= n_times; i++){

      relative_error[i][0] = abs((analytical_values[0]-expectation_values[i][0])/analytical_values[0]);   //Stores relative error in E
      relative_error[i][1] = abs((analytical_values[1]-expectation_values[i][1])/analytical_values[1]);   //Stores relative error in E_squared
      relative_error[i][2] = abs((analytical_values[2]-expectation_values[i][2])/analytical_values[2]);   //Stores relative error in Mabs
      relative_error[i][3] = abs((analytical_values[3]-expectation_values[i][3])/analytical_values[3]);   //Stores relative error in Mabs_squared
      relative_error[i][4] = abs((analytical_values[5]-expectation_values[i][5])/analytical_values[5]);   //Stores relative error in M_squared
      relative_error[i][5] = abs((analytical_values[6]-expectation_values[i][6])/analytical_values[6]);   //Stores relative error in Heat Capacity
      relative_error[i][6] = abs((analytical_values[7]-expectation_values[i][7])/analytical_values[7]);   //Stores relative error in Magnetic susceptibility

    ofile << setprecision(9) << time[i] << " " << setprecision(9) << expectation_values[i][0] << " " << setprecision(9) << expectation_values[i][1] << " " << setprecision(9) << expectation_values[i][2]
    << " " << setprecision(9) << expectation_values[i][3]<< " " << setprecision(9) << expectation_values[i][4] << " " << setprecision(9) << expectation_values[i][5] << " " << setprecision(9) << expectation_values[i][6]
    << " " << setprecision(9) << expectation_values[i][7] << endl ;


    ofile2 << setprecision(9) << time[i] << " " << setprecision(9) << relative_error[i][0] << " " << setprecision(9) << relative_error[i][1] << " " << setprecision(9) << relative_error[i][2]
    << " " << setprecision(9) << relative_error[i][3]<< " " << setprecision(9) << relative_error[i][4] << " " << setprecision(9) << relative_error[i][5] << " " << setprecision(9) << relative_error[i][6] << endl;
    }

    ofile.close();
    ofile2.close();
}

return 0;
}

void initialize(int dimensions, int **spin_matrix, double& E, double& M, string initialize){
  if (initialize == "ordered"){
    // setup spin matrix and intial magnetization
    for(int i =0; i <dimensions; i++) {
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

  }
  // setup initial energy
  for(int i = 0; i < dimensions; i++) {
    for (int j= 0; j < dimensions; j++){
      E -=  (double) spin_matrix[i][j] * (spin_matrix[periodic(i,dimensions,-1)][j] + spin_matrix[i][periodic(j,dimensions,-1)]);
      }
    }
  }

void analytical_values_2x2Lattice(double* analytical_values, double T){

  double Z_a,E_a,M_a,E_squared_a,M_squared_a, Mabs_a, Mabs_squared_a,beta;

  beta = 1/T;

  //Calculates the analytical expectation values
  Z_a =  4*(3 + cosh(8*beta));
  E_a = -32*sinh(8*beta)/Z_a;
  E_squared_a = 256*cosh(8*beta)/Z_a;
  M_a = 0;
  M_squared_a = 32*(exp(8*beta) + 1)*beta/Z_a ;
  Mabs_a =  8*(exp(8*beta) + 2)/Z_a;
  Mabs_squared_a = (32*exp(8*beta) + 4)/Z_a;


  //Stores the analytical expectation values
  analytical_values[0] = E_a;
  analytical_values[1] = E_squared_a;
  analytical_values[2] = Mabs_a;
  analytical_values[3] = Mabs_squared_a;
  analytical_values[4] = M_a;
  analytical_values[5] = M_squared_a;
  analytical_values[6] = (analytical_values[1]-analytical_values[0]*analytical_values[0])*beta*beta;    //Heat Capacity
  analytical_values[7] = (analytical_values[5])*beta;                                                   //Magnetic susceptibility


}

void Monte_Carlo_Metropolis_time(int MC, int n, int **spin_matrix, int J, double& E, double& M, double* boltzmann_distribution, double* energy, double* magnetization, double* time, double* acceptance, double beta){


  //SPØR OM DET HER KAN SENDES INN SÅ DET IKKE MÅ LAGES HVER ITERASJON!!!!!
  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntegerGenerator(0,n-1);        //Sets up the uniform distribution for x in [0,n-1]
  uniform_real_distribution<double> RandomNumberGenerator(0,1);       //Sets up the uniform distribution for x in [0,1]

  int x_flip, y_flip, dE, dM, n_spins, i, accept;
  double E_sum, M_sum, Mabs_sum, Mabs_sum_squared, E_squared, M_squared;

  n_spins = n*n;

  E_sum = 0.0;
  M_sum = 0.0;
  Mabs_sum = 0.;
  Mabs_sum_squared = 0.;
  E_squared = E*E;
  M_squared = M*M;
  accept = 0;


  //Running over Monte Carlo samples
  for (int k = 0; k < MC; k++){


    //for (int j = 0; j < n_spins; j++){
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
        accept += 1;
      }
      else if(RandomNumberGenerator(gen) < boltzmann_distribution[dE + 8]){
        spin_matrix[x_flip][y_flip] *= (-1);     //Accepting flip
        dM = 2*spin_matrix[x_flip][y_flip];
        accept += 1;
      }
      else{
        dE = 0;                             //Rejecting the flip
        dM = 0;
      }

      E += dE;
      M += dM;

    //}

    E_sum += (double) E;
    M_sum += (double) M;
    Mabs_sum += (double) abs(M);

    E_squared += (double) E*E;
    M_squared += (double) M*M;
    Mabs_sum_squared += (double) abs(M)*abs(M);

    if (k % n_spins == 0){
      i = k / n_spins;
      energy[i] = E_sum/ (( k * n_spins);
      magnetization[i] = M_sum / ( k * n_spins);
      time[i] = i;
      acceptance[i] = accept;
    }
  }



  E_squared /= (double) MC;
  M_squared /= (double) MC;
  Mabs_sum /= (double) MC;
  Mabs_sum_squared /= (double) MC;


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

void Monte_Carlo_Metropolis_2x2(int MC, int n, int **spin_matrix, int J, double& E, double& M, double* boltzmann_distribution, double* time, double* acceptance, double beta,double** expectation_values){


  //SPØR OM DET HER KAN SENDES INN SÅ DET IKKE MÅ LAGES HVER ITERASJON!!!!!
  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntegerGenerator(0,n-1);        //Sets up the uniform distribution for x in [0,n-1]
  uniform_real_distribution<double> RandomNumberGenerator(0,1);       //Sets up the uniform distribution for x in [0,1]

  int x_flip, y_flip, dE, dM, n_spins, i, accept, n_times;
  double E_sum, M_sum, Mabs_sum, Mabs_sum_squared, E_sum_squared, M_sum_squared;
  double* Mabs, Mabs_squared, E_squared, M_squared, energy, magnetization;

  n_spins = n*n;
  n_times = MC_samples/n_spins;

  E_sum = 0.0;
  M_sum = 0.0;
  Mabs_sum = 0.;
  Mabs_sum_squared = 0.;
  E_sum_squared = E*E;
  M_sum_squared = M*M;
  accept = 0;

  energy = new double[n_times];
  magnetization = new double[n_times];
  Mabs = new double[n_times];
  Mabs_squared = new double[n_times];
  E_squared = new double[n_times];
  M_squared = new double[n_times];


  //Running over Monte Carlo samples
  for (int k = 0; k < MC; k++){


    //for (int j = 0; j < n_spins; j++){
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
        accept += 1;
      }
      else if(RandomNumberGenerator(gen) < boltzmann_distribution[dE + 8]){
        spin_matrix[x_flip][y_flip] *= (-1);     //Accepting flip
        dM = 2*spin_matrix[x_flip][y_flip];
        accept += 1;
      }
      else{
        dE = 0;                             //Rejecting the flip
        dM = 0;
      }

      E += dE;
      M += dM;

    //}

    E_sum += (double) E;
    M_sum += (double) M;
    Mabs_sum += (double) abs(M);

    E_sum_squared += (double) E*E;
    M_sum_squared += (double) M*M;
    Mabs_sum_squared += (double) abs(M)*abs(M);

    if (k % n_spins == 0){
      i = k / n_spins;
      energy[i] = E_sum/ (k * n_spins);
      magnetization[i] = M_sum / ( k * n_spins);
      Mabs[i] = Mabs_sum/ ( k * n_spins);
      Mabs_squared[i] = Mabs_sum_squared/ ( k * n_spins);
      E_squared[i] = E_sum_squared/ ( k * n_spins);
      M_squared[i] = M_sum_squared/ ( k * n_spins);
      time[i] = i;
      acceptance[i] = accept;
    }



  }
  expectation_values[0] = energy;
  expectation_values[1] = E_squared;
  expectation_values[2] = Mabs;
  expectation_values[3] = Mabs_squared;
  expectation_values[4] = magnetization;
  expectation_values[5] = M_squared;
  expectation_values[6] = (expectation_values[1]-expectation_values[0]*expectation_values[0])*beta*beta;             //Stores the computed expectation value for heat capacity
  expectation_values[7] = (expectation_values[5]-(expectation_values[4]*expectation_values[4]))*beta;      //Stores the computed expectation value for susceptibility

}
