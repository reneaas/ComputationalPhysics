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
using namespace std;

ofstream ofile;


inline int periodic(int coordinate, int dimensions, int step) {
  return (coordinate+dimensions+step) % (dimensions);
}


void initialize(int , int **, double& , double& , string);
void Monte_Carlo_Metropolis_time(int, int, int, int**, int, double&, double&, double&, double&, double&, double& , double* , double , double& , mt19937_64 , uniform_int_distribution<int> , uniform_real_distribution<double> );

int main(int nargs, char* args[]){
  string outfilename, local_outfilename;
  double boltzmann_factors[17];
  int L, MC_samples, J;
  int **spin_matrix, my_rank, numprocs, **local_spin_matrix;
  double E_initial, M_initial;                      //Stores initial energy and magnetization of system.
  int n_spins;                                      //Total number of spins.
  double beta;
  int local_n;                                      //Number of local temperature steps.
  double h, local_h;                                         //step size of temperature.
  double local_T0, local_T1, T_start, T_final;
  int number_of_temperatures, N;
  double temp;
  double E_local, M_local, Cv_local, chi_local, variance_local, time_start, time_end, timeused;


  //Read from command line
  L = atoi(args[1]);                                    //Dimension of spin matrix
  MC_samples = atoi(args[2]);                           //Number of Monte Carlo samples
  N = atoi(args[3]);                                    //Burn-in period.
  number_of_temperatures = 80;
  //initialize matrix
  spin_matrix = new int*[L];
  for (int i = 0; i < L; i++){
    spin_matrix[i] = new int[L];
  }

  E_initial = 0;
  M_initial = 0;
  n_spins = L*L;
  J = 1;                                                      //Coupling constant
  initialize(L, spin_matrix, E_initial, M_initial, "ordered");

  // MPI initializations
  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  ofstream local_ofile;
  local_n = number_of_temperatures/numprocs;
  T_start = 2.1;
  T_final = 2.4;
  E_local = 0.; M_local = 0.; Cv_local = 0.; chi_local = 0.;
  h = (T_final - T_start)/((double) number_of_temperatures);
  local_T0 = T_start + (double) my_rank*local_n*h;
  local_T1 = local_T0 + (double) local_n*h;
  local_h = (local_T1 - local_T0)/((double) local_n);
  local_outfilename = "observables_my_rank_" + to_string(my_rank) + "_L_" + to_string(L) + ".txt";
  //cout << "h local = " << local_h << endl;
  local_spin_matrix = spin_matrix;

  random_device rd;
  mt19937_64 gen(rd() + my_rank);
  uniform_int_distribution<int> RandomIntegerGenerator(0,L-1);        //Sets up the uniform distribution for x in [0,n-1]
  uniform_real_distribution<double> RandomNumberGenerator(0,1);       //Sets up the uniform distribution for x in [0,1]

  time_start = MPI_Wtime();

  local_ofile.open(local_outfilename);
  for (int i = 0; i < local_n; i++){
    if (my_rank == 0){
      cout << "iteration = " << i << " of " << number_of_temperatures/numprocs << endl;
    }
    temp = local_T0 + (double) i*local_h;
    //Compute Boltzmann factors.
    beta = 1/(temp);                //k_B = 1
    for (int j = -8; j < 9; j += 4){
      boltzmann_factors[j + 8] = exp(-beta*j);
    }
    local_spin_matrix = spin_matrix;
    E_local = E_initial;
    M_local = M_initial;

    Monte_Carlo_Metropolis_time(MC_samples, L, N, local_spin_matrix, J, E_initial, E_local, M_initial, M_local, chi_local, Cv_local, boltzmann_factors, beta, variance_local, gen, RandomIntegerGenerator, RandomNumberGenerator);
    //MPI_Barrier (MPI_COMM_WORLD);

    //cout << "my_rank = " << my_rank << " " << "E = " << E_local << endl;

    local_ofile << temp << " " << E_local << " " << M_local << " " << chi_local << " " << Cv_local << endl;
    E_local = 0.;
    M_local = 0.;
    chi_local = 0.;
    Cv_local = 0.;

  }

  local_ofile.close();

  time_end = MPI_Wtime();
  timeused = time_end - time_start;
  if (my_rank == 0){
    cout << "Time used = " << timeused << " seconds" <<  endl;
  }

  MPI_Finalize();

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

void Monte_Carlo_Metropolis_time(int MC, int n, int N, int **spin_matrix, int J, double& E, double& E_exp, double& M, double& M_exp, double& chi, double& Cv, double* boltzmann_factors, double beta, double& variance, mt19937_64 gen, uniform_int_distribution<int> RandomIntegerGenerator, uniform_real_distribution<double> RandomNumberGenerator){


    int x_flip, y_flip, dE, dM, n_spins, i;
    double E_sum, M_sum, Mabs_sum, Mabs_sum_squared, E_squared, M_squared;

    n_spins = n*n;
    E_sum = 0.0;
    M_sum = 0.0;
    Mabs_sum = 0.;
    Mabs_sum_squared = 0.;
    E_squared = 0.;
    M_squared = 0.;


    //Running over Monte Carlo samples
    for (int k = 1; k <= MC; k++){


      x_flip = RandomIntegerGenerator(gen);
      y_flip = RandomIntegerGenerator(gen);      //Randomized indices to matrix element that will be flipped


      dE = (int) 2*J*spin_matrix[x_flip][y_flip] * (spin_matrix[periodic(x_flip,n,1)][y_flip] + spin_matrix[periodic(x_flip, n,-1)][y_flip]
                                                    + spin_matrix[x_flip][periodic(y_flip, n,1)] + spin_matrix[x_flip][periodic(y_flip, n,-1)]);

      //Metropolis algorithm
      if(dE < 0){
        spin_matrix[x_flip][y_flip] *= (-1);   //Accepting the flip
        dM = 2*spin_matrix[x_flip][y_flip];
      }
      else if(RandomNumberGenerator(gen) < boltzmann_factors[dE + 8]){
        spin_matrix[x_flip][y_flip] *= (-1);     //Accepting flip
        dM = 2*spin_matrix[x_flip][y_flip];
      }
      else{
        dE = 0;                             //Rejecting the flip
        dM = 0;
      }

      E += (double) dE;
      M += (double) dM;


      if (k > N-1){
        E_sum += (double) E;
        M_sum += (double) M;
        Mabs_sum += (double) abs(M);

        E_squared += (double) E*E;
        M_squared += (double) M*M;
        Mabs_sum_squared += (double) abs(M)*abs(M);
      }

    }

    E_sum /= (double) (MC-N);
    E_squared /= (double) (MC-N);
    variance = E_squared - E_sum*E_sum;
    Cv = variance*beta*beta/((double) n_spins);

    M_sum /= (double) (MC-N);
    M_squared /= (double) (MC-N);

    Mabs_sum /= (double) (MC-N);
    Mabs_sum_squared /= (double) (MC-N);
    chi = (Mabs_sum_squared - Mabs_sum*Mabs_sum)*beta/((double) n_spins);


    E_exp = E_sum/((double) n_spins);
    M_exp = Mabs_sum/((double) n_spins);
    //E = E_sum; M = Mabs_sum;

}
