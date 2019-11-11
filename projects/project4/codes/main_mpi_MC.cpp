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
void Monte_Carlo_Metropolis_time( int, int, int, int**, int, double&, double&, double&, double&, double&, double& , double* , double , mt19937_64 , uniform_int_distribution<int> , uniform_real_distribution<double> );

int main(int nargs, char* args[]){
  string outfilename;
  double boltzmann_factors[17];
  int L, J;
  int MC, MC_local;
  int  n_spins, N_local;
  int **spin_matrix, my_rank, numprocs, **local_spin_matrix;
  double E_initial, M_initial;                      //Stores initial energy and magnetization of system.
  double beta, chi, Cv;
  double E_sum, M_sum, Esq_sum, Msq_sum;
  double E_local, M_local, Esum_local, Msum_local, Esq_local, Msq_local;
  double T_init, T_final, h;
  double time_start, time_end, timeused, N, number_of_temperatures, temp;


  //Read from command line
  L = atoi(args[1]);                                    //Dimension of spin matrix
  MC = atol(args[2]);                           //Number of Monte Carlo samples
  N = atoi(args[3]);                                    //Burn-in period.
  outfilename = "observables_L_" + to_string(L) + ".txt";

  //initialize matrix
  spin_matrix = new int*[L];
  for (int i = 0; i < L; i++){
    spin_matrix[i] = new int[L];
  }
  E_initial = 0;
  M_initial = 0;                                                 //Coupling constant
  initialize(L, spin_matrix, E_initial, M_initial, "random");



  // MPI initializations
  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  E_sum = 0;
  M_sum = 0;
  Esq_sum = 0;
  Msq_sum = 0;
  n_spins = L*L;
  J = 1;
  T_init = 2.0;
  T_final = 2.4;
  number_of_temperatures = 400;
  h = (T_final - T_init)/number_of_temperatures;
  MC_local = MC/numprocs;
  N_local = N/numprocs;

  if (my_rank == 0){
    cout << "local MC = " << MC_local << endl;
    cout << "Burn in period = " << N << endl;
  }

  random_device rd;
  mt19937_64 gen(rd() + my_rank);
  uniform_int_distribution<int> RandomIntegerGenerator(0,L-1);        //Sets up the uniform distribution for x in [0,n-1]
  uniform_real_distribution<double> RandomNumberGenerator(0,1);       //Sets up the uniform distribution for x in [0,1]

  time_start = MPI_Wtime();
  if (my_rank == 0) ofile.open(outfilename);

  for (int i = 0; i < number_of_temperatures; i++){
    temp = T_init + (double) i*h;
    if (my_rank == 0) cout << "T = " << temp << endl;
    //Compute Boltzmann factors.
    beta = 1/(temp);                //k_B = 1
    for (int j = -8; j < 9; j += 4){
      boltzmann_factors[j + 8] = exp(-beta*j);
    }

    local_spin_matrix = spin_matrix;
    E_local = E_initial;
    M_local = M_initial;


    Monte_Carlo_Metropolis_time(MC_local, L, N_local, local_spin_matrix, J, E_initial, Esum_local, M_initial, Msum_local, Esq_local, Msq_local, boltzmann_factors, beta, gen, RandomIntegerGenerator, RandomNumberGenerator);

    MPI_Reduce(&Esum_local, &E_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Msum_local, &M_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Esq_local, &Esq_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Msq_local, &Msq_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    E_sum /= (double) (MC-N); M_sum /= (double) (MC-N); Esq_sum /= (double) (MC-N); Msq_sum /= (double) (MC-N);

    chi = (Msq_sum - M_sum*M_sum)*beta / n_spins ;
    Cv = (Esq_sum - E_sum*E_sum)*beta*beta / n_spins;

    E_sum /= n_spins; M_sum /= n_spins; Esq_sum /= n_spins; Msq_sum /= n_spins;

    ofile << temp << " " << E_sum << " " << M_sum << " " << chi << " " << Cv << endl;

    Esum_local = 0.; E_sum = 0; Esq_local = 0; Esq_sum = 0;
    Msum_local = 0.; M_sum = 0; Msq_local = 0; Msq_sum = 0;

  }

  if (my_rank == 0) ofile.close();

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

void Monte_Carlo_Metropolis_time( int MC, int n, int N, int **spin_matrix, int J, double& E, double& E_exp, double& M, double& M_exp, double& Esq_exp, double& Msq_exp, double* boltzmann_factors, double beta, mt19937_64 gen, uniform_int_distribution<int> RandomIntegerGenerator, uniform_real_distribution<double> RandomNumberGenerator){


    int x_flip, y_flip, dE, dM;
    double E_sum, M_sum, Mabs_sum, Mabs_sum_squared, E_squared, M_squared;

    E_sum = 0.0;
    M_sum = 0.0;
    Mabs_sum = 0.;
    Mabs_sum_squared = 0.;
    E_squared = 0.;
    M_squared = 0.;


    //Running over Monte Carlo samples
    for ( int k = 1; k <= MC; k++){
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

    E_exp = E_sum;
    M_exp = Mabs_sum;
    Esq_exp = E_squared;
    Msq_exp = Mabs_sum_squared;
}
