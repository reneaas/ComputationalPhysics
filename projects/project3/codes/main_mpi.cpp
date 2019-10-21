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
#define EPS 3.0e-14
#define MAXIT 10

using namespace std;

ofstream ofile;

double CoulombRepulsionMC_ImportanceSampling(double*);
double CoulombRepulsionMC(double*);
double laguerre_integrate_func_6(double, double, double, double, double, double);
double gammln(double);
void gauss_laguerre(double*, double*, int, double);


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

  if (integration_method == "3"){
    int n, local_n, numprocs, my_rank;
    double total_integral, local_integral;
    double a, local_a, b, local_b, h, alpha;
    double time_start, time_end, total_time;

    double exact, rel_error;
    double *w1, *w2, *w3, *w4, *w5, *w6, *r1, *r2, *theta1, *theta2, *phi1, *phi2;
    int local_n2;

    a = 0;
    b = M_PI;
    n = 20;

    total_integral = 0.;
    local_integral = 0.;
    alpha = 0;
    h = (b-a)/((double) n);
    r1 = new double[n+1];
    r2 = new double[n+1];
    w1 = new double[n+1];
    w2 = new double[n+1];


    gauss_laguerre(r1, w1, n, alpha);
    gauss_laguerre(r2, w2, n, alpha);


    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    local_n = n/numprocs;
    w3 = new double[local_n];
    w4 = new double[local_n];
    w5 = new double[local_n];
    w6 = new double[local_n];


    theta1 = new double[local_n];
    theta2 = new double[local_n];

    phi1 = new double[local_n];
    phi2 = new double[local_n];

    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;

    gauleg(local_a, local_b, theta1, w3, local_n);
    gauleg(local_a, local_b, theta2,w4, local_n);
    gauleg(local_a, 2*local_b, phi1, w5, local_n);
    gauleg(local_a, 2*local_b, phi2, w6, local_n);

    local_n2 = local_n*my_rank;



    time_start = MPI_Wtime();


    for (int i = 1 ; i <= n; i++){
      if (my_rank == 0){
        cout << "iteration = " << i << endl;
      }
      for (int j = 1; j <= n; j++){
        for (int k = 0; k < local_n; k++){
          for (int l = 0; l < local_n; l++){
            for (int p = 0; p < local_n; p++){
              for (int r = 0; r < local_n; r++){
                local_integral += w1[i]*w2[j]*w3[k]*w4[l]*w5[p]*w6[r]*laguerre_integrate_func_6(r1[i], r2[j], theta1[k], theta2[l], phi1[p], phi2[r]);
              }
            }
          }
        }
      }
    }




    MPI_Reduce(&local_integral, &total_integral, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    time_end = MPI_Wtime();

    total_time = time_end - time_start;
    exact = 5*pow(M_PI, 2)/(16*16);
    rel_error = abs(total_integral - exact)/exact;

    if (my_rank == 0){
    cout << "Exact = " << " " << exact << endl;
    cout << "Calculated = " << " " <<total_integral << endl;
    cout << "Total time =  " << total_time << endl;
    }
    /*
    if (my_rank == 0){
    ofile.open(outfilename);
    ofile << n << " " << setprecision(9) << total_integral << " "<< setprecision(9) << rel_error << " " <<setprecision(9) <<total_time <<endl;
    ofile.close();
    */

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

double laguerre_integrate_func_6(double u1, double u2, double theta1, double theta2, double phi1, double phi2){
  if (u1*u1 + u2*u2 - 2*u1*u2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1 - phi2)) <  1e-8){
    return 0;
  }

  return exp(-3*(u1+u2))*(sin(theta2)*sin(theta1)*u1*u1*u2*u2/sqrt(u1*u1 + u2*u2 - 2*u1*u2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1 - phi2))));

}

void gauss_laguerre(double *x, double *w, int n, double alf){
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}
// end function gaulag

double gammln( double xx){
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln
#undef EPS
#undef MAXIT
