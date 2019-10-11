#include <iostream>
#include "lib.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <string>
#include "time.h"
#include <random>
#include <fstream>
#define EPS 3.0e-14
#define MAXIT 10


using namespace std;

ofstream ofile;

double CoulombRepulsion(double, double, double, double, double, double);
double CoulombRepulsion_spherical(double, double, double);
double F(double, double);
double radial_probability_density(double,double);
double test_func(double, double, double);
double gammln(double);
void gauss_laguerre(double*, double*, int, double);
double laguerre_integrate_func_3(double, double, double);
double laguerre_integrate_func_6(double, double, double, double, double, double);

int main(int nargs, char* args[]){
  string integration_method;
  string outfilename;
  if (nargs == 1){

    cout << "Specify integration method, choose from: " << endl;
    cout << "------------------------------------------------------" << endl;
    cout << "Gauss Legendre method                  --> type 1 " << endl;
    cout << "Gauss Laguerre method                  --> type 2 " << endl;
    cout << "Brute force Monte Carlo                --> type 3 " << endl;
    cout << "Monte Carlo with importance sampling   --> type 4 " << endl;
    cout << "------------------------------------------------------" << endl;
    cin >> integration_method;
  }
  else{
    integration_method = "montecarlo_benchmarking";
    outfilename = string(args[1]);
  }

  if (integration_method == "1"){
    int n;
    double a,b;
    cout << "Specify number of integration points: " << endl;
    cin >> n;
    cout << "Specify integration limits [a,b] : " << endl;
    cin >> a >> b;

    double *x1, *y1, *z1, *x2, *y2, *z2;
    double *w1, *w2, *w3, *w4, *w5, *w6;
    x1 = new double[n];
    y1 = new double[n];
    z1 = new double[n];
    x2 = new double[n];
    y2 = new double[n];
    z2 = new double[n];

    w1 = new double[n];
    w2 = new double[n];
    w3 = new double[n];
    w4 = new double[n];
    w5 = new double[n];
    w6 = new double[n];

    //Generate mesh points and weights using the Gauss-Legendre method found in Numerical Recipies.
    gauleg(a,b,x1,w1,n);
    gauleg(a,b,y1,w2,n);
    gauleg(a,b,z1,w3,n);
    gauleg(a,b,x2,w4,n);
    gauleg(a,b,y2,w5,n);
    gauleg(a,b,z2,w6,n);

    //Time the program:
    clock_t start, finish;
    start = clock();
    double integral_gauss_legendre = 0;
    //Solves the integral using the weights and mesh points obtained with the library function gauleg.
    for (int i = 0; i < n; i++){
      for (int j = 0; j < n; j++){
        for (int k = 0; k < n; k++){
          for (int l = 0; l < n; l++){
            for (int p = 0; p < n; p++){
              for (int r = 0; r < n; r++){
                integral_gauss_legendre += w1[i]*w2[j]*w3[k]*w4[l]*w5[p]*w6[r]*CoulombRepulsion(x1[i], y1[j], z1[k], x2[l], y2[p], z2[r]);
              }
            }
          }
        }
      }
    }
    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC);
    cout << "Integral = " << integral_gauss_legendre << endl;
    cout << "Analytical value = " << 5*pow(M_PI, 2)/(16*16) << endl;
    cout << "timeused = " << timeused << endl;
  }

  if (integration_method == "2"){
    int n, dimensions;
    double a,b, alpha;
    cout << "Specify number of integration points: " << endl;
    cin >> n;

    a = 0;
    b = M_PI;
    double integral_gauss_laguerre = 0;


    cout<<"Choose dimension for integral:"<<endl;
    cout << "------------------------------------------------------" << endl;
    cout<<"For 3 dimensions                         --> type 3 "<<endl;
    cout<<"For 6 dimensions                         --> type 6 "<<endl;
    cout << "------------------------------------------------------" << endl;
    cin >> dimensions;
    cout << "------------------------------------------------------" << endl;

    if(dimensions == 3){
      double *w1, *w2, *w3, *r1, *r2, *theta2;

    w1 = new double[n];
    w2 = new double[n];
    w3 = new double[n];
    r1 = new double[n+1];
    r2 = new double[n+1];
    theta2 = new double[n];

    alpha = 2.0;

    gauss_laguerre(r1, w1, n, alpha);
    gauss_laguerre(r2, w2, n, alpha);
    gauleg(a, b, theta2, w3, n);


    for(int i = 1; i<(n+1); i++){
      for(int j = 1; j<(n+1); j++){
        for(int k = 0; k<n; k++){
          integral_gauss_laguerre += w1[i]*w2[j]*w3[k]*laguerre_integrate_func_3(r1[i], r2[j], theta2[k]);
        }
      }
    }

    integral_gauss_laguerre *= (M_PI*M_PI)/128;
    }

    if(dimensions == 6){
      double *w1, *w2, *w3, *w4, *w5, *w6, *r1, *r2, *theta1, *theta2, *phi1, *phi2;

      w1 = new double[n+1];
      w2 = new double[n+1];
      w3 = new double[n];
      w4 = new double[n];
      w5 = new double[n];
      w6 = new double[n];

      r1 = new double[n+1];
      r2 = new double[n+1];

      theta1 = new double[n];
      theta2 = new double[n];

      phi1 = new double[n];
      phi2 = new double[n];

      alpha = 0;

      gauss_laguerre(r1, w1, n, alpha);
      gauss_laguerre(r2, w2, n, alpha);
      gauleg(a, b, theta1, w3, n);
      gauleg(a, b, theta2,w4, n);
      gauleg(a, 2*b, phi1, w5, n);
      gauleg(a, 2*b, phi2, w6, n);



      for (int i = 1; i <= n; i++){
        for (int j = 1; j <= n; j++){
          for (int k = 0; k < n; k++){
            for (int l = 0; l < n; l++){
              for (int p = 0; p < n; p++){
                for (int r = 0; r < n; r++){
                  integral_gauss_laguerre += w1[i]*w2[j]*w3[k]*w4[l]*w5[p]*w6[r]*laguerre_integrate_func_6(r1[i], r2[j], theta1[k], theta2[l], phi1[p], phi2[r]);
                }
              }
            }
          }
        }
      }


    }


    cout << "Integral = " << integral_gauss_laguerre << endl;
    cout << "Analytical value = " << 5*pow(M_PI, 2)/(16*16) << endl;


  }



  if (integration_method == "3"){
    int n;
    double a,b;                   // The integration interval [a,b].
    int N;                //number of Monte Carlo samples
    cout << "Read in the number of integration points" << endl;
    cin >> n;
    cout << "Read in the integration limits [a,b] " << endl;
    cin >> a >> b;
    cout << "Specify number of monte carlo samples: " << endl;
    cin >> N;
    //Initialize the seed and call the Mersienne algorithm
    int d = 6;          //d-dimensional integral.
    random_device rd;
    mt19937_64 gen(rd());
    double x1, x2, x3, x4, x5, x6;
    //Sets up the uniform distribution for x in [0,1]
    uniform_real_distribution<double> RandomNumberGenerator(0,1);
    //Creates the variables
    double MC_integral;
    double *MC_integrals;
    double *x = new double[d];
    MC_integrals = new double[N];
    for (int i = 0; i < N; i++){
      cout << "Computing for sample = " << i << endl;
      for (int j = 0; j < n; j++){
        x[0] =  RandomNumberGenerator(gen);
        x[0] = a + (b-a)*x[0];
        for (int k = 0; k < n; k++){
          x[1] = RandomNumberGenerator(gen);
          x[1] = a + (b-a)*x[1];
          for (int l = 0; l < n; l++){
            x[2] = RandomNumberGenerator(gen);
            x[2] = a + (b-a)*x[2];
            for (int m = 0; m < n; m++){
              x[3] = RandomNumberGenerator(gen);
              x[3] = a + (b-a)*x[3];
              for (int p = 0; p < n; p++){
                x[4] = RandomNumberGenerator(gen);
                x[4] = a + (b-a)*x[4];
                for (int r = 0; r < n; r++){
                  x[5] = RandomNumberGenerator(gen);
                  x[5] = a + (b-a)*x[5];
                  MC_integrals[i] += CoulombRepulsion(x[0],x[1],x[2],x[3],x[4],x[5]);
                }
              }
            }
          }
        }
      }
      MC_integral += MC_integrals[i]/ (pow((double) n,d));
    }
    MC_integral /= (double) N;
    MC_integral *= pow((double) (b-a), d);      //Compensates for the change of variables xi = a + (b-a)*mu.

    cout << "Computed integral = " << MC_integral << endl;
    cout << "Analytical value = " << 5*pow(M_PI,2)/(16*16) << endl;
  }

  if (integration_method == "4"){
    int n;
    double max_radial_distance;
    cout << "Read in the number of integration points" << endl;
    cin >> n;
    cout << "Read in maximum radial distance r " << endl;
    cin >> max_radial_distance;
    int d = 6;          //d-dimensional integral.
    double r1, r2, theta2;
    double *a, *b;
    a = new double[d];
    b = new double[d];
    double alpha = 4;

    //r1 endpoints
    a[0] = 0;
    b[0] = max_radial_distance;
    //r2 endpoints
    a[1] = 0;
    b[1] = max_radial_distance;
    //theta1 endpoints
    a[2] = 0;
    b[2] = M_PI;
    //theta2 endpoints
    a[3] = 0;
    b[3] = M_PI;
    //phi1 endpoints
    a[4] = 0;
    b[4] = 2*M_PI;
    //phi2 endpoints
    a[5] = 0;
    b[5] = 2*M_PI;

    //Initialize the seed and call the Mersienne algorithm
    random_device rd;
    mt19937_64 gen(rd());
    //Sets up the uniform distribution for x in [0,1]
    uniform_real_distribution<double> RandomNumberGenerator(0,1);
    //Creates the variables
    double MC_integral = 0;
    double *MC_integrals;
    int N;               //number of Monte Carlo samples
    cout << "Specify the number of Monte Carlo samples: " << endl;
    cin >> N;
    MC_integrals = new double[N];

    //Benchmark code
    clock_t start, finish;
    start  = clock();
    //Main algorithm
    for (int i = 0; i < N; i++){
      cout << "Computing for sample = " << i << endl;
      for (int j = 0; j < n; j++){
        r1 = RandomNumberGenerator(gen);
        r1 = -log(1-r1)/alpha;
        for (int k = 0; k < n; k++){
          r2 = RandomNumberGenerator(gen);
          r2 = -log(1-r2)/alpha;
          for (int l = 0; l < n; l++){
            theta2 = RandomNumberGenerator(gen);
            theta2 = a[3] + (b[3]-a[3])*theta2;
            MC_integrals[i] += CoulombRepulsion_spherical(r1,r2,theta2)/radial_probability_density(r1,r2);
          }
        }
      }
      MC_integral += MC_integrals[i] / (pow( (double) n, (double) 3));
    }
    MC_integral /= (double) N;
    MC_integral *= 8*M_PI*M_PI;                               //Multiplying by factors due to integration with respect to phi1, phi2 and theta1. Integrand not explicitly dependent on them.
    finish = clock();
    double timeused = (double) (finish-start)/(CLOCKS_PER_SEC);
    MC_integral *= M_PI;                            //Since we're using a uniform distribution for theta2 only, we only need to multiply by pi.


    cout << "Computed integral = " << MC_integral << endl;
    cout << "Analytical value = " << 5*pow(M_PI,2)/(16*16) << endl;
    cout << "time used = " << timeused << endl;
  }

  if (integration_method == "montecarlo_benchmarking"){
    int n = 10;
    double max_radial_distance = 10;
    int d = 6;          //d-dimensional integral.
    double r1, r2, theta2;
    double *a, *b;
    a = new double[d];
    b = new double[d];
    double alpha = 4;

    //r1 endpoints
    a[0] = 0;
    b[0] = max_radial_distance;
    //r2 endpoints
    a[1] = 0;
    b[1] = max_radial_distance;
    //theta1 endpoints
    a[2] = 0;
    b[2] = M_PI;
    //theta2 endpoints
    a[3] = 0;
    b[3] = M_PI;
    //phi1 endpoints
    a[4] = 0;
    b[4] = 2*M_PI;
    //phi2 endpoints
    a[5] = 0;
    b[5] = 2*M_PI;

    //Initialize the seed and call the Mersienne algorithm
    random_device rd;
    mt19937_64 gen(rd());
    //Sets up the uniform distribution for x in [0,1]
    uniform_real_distribution<double> RandomNumberGenerator(0,1);
    //Creates the variables
    double MC_integral = 0;
    double *MC_integrals;
    int N;               //number of Monte Carlo samples
    N = atoi(args[2]);
    MC_integrals = new double[N];

    //Benchmark code
    clock_t start, finish;
    start  = clock();
    //Main algorithm
    for (int i = 0; i < N; i++){
      //cout << "Computing for sample = " << i << endl;
      for (int j = 0; j < n; j++){
        r1 = RandomNumberGenerator(gen);
        r1 = -log(1-r1)/alpha;
        for (int k = 0; k < n; k++){
          r2 = RandomNumberGenerator(gen);
          r2 = -log(1-r2)/alpha;
          for (int l = 0; l < n; l++){
            theta2 = RandomNumberGenerator(gen);
            theta2 = a[3] + (b[3]-a[3])*theta2;
            MC_integrals[i] += CoulombRepulsion_spherical(r1,r2,theta2)/radial_probability_density(r1,r2);
          }
        }
      }
      MC_integral += MC_integrals[i] / (pow( (double) n, (double) 3));
    }
    MC_integral /= (double) N;
    MC_integral *= 8*M_PI*M_PI;                               //Multiplying by factors due to integration with respect to phi1, phi2 and theta1. Integrand not explicitly dependent on them.
    finish = clock();
    double timeused = (double) (finish-start)/(CLOCKS_PER_SEC);
    MC_integral *= M_PI;                            //Since we're using a uniform distribution for theta2 only, we only need to multiply by pi.

    /*
    cout << "Computed integral = " << MC_integral << endl;
    cout << "Analytical value = " << 5*pow(M_PI,2)/(16*16) << endl;
    cout << "time used = " << timeused << endl;
    */

    ofile.open(outfilename);
    ofile << N << " " << timeused << " " << MC_integral << endl;
    ofile.close();

  }

  return 0;
}

double CoulombRepulsion(double x1, double y1, double z1, double x2, double y2, double z2){
  /*
  The function which we integrate.
  */
  if (sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)) == 0){
    return 0;
  }
  return exp( -4*( sqrt(x1*x1 + y1*y1 + z1*z1) + sqrt(x2*x2 + y2*y2 + z2*z2) ))/sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
}

double CoulombRepulsion_spherical(double r1, double r2, double theta2){
  /*
  CoulombRepulsion in spherical coordinates. Includes the Jacobi determinant r1^2r2^2 sin(theta2).
  The other coordinates are just multiplied as constants since the integral is not explicitly dependent upon phi1, phi2 or theta1.
  */
  if (sqrt(r1*r1 + r2*r2 - r1*r2*cos(theta2)) == 0){
    return 0;
  }
  return r1*r1*r2*r2*sin(theta2)*exp(-4*(r1+r2)) / ( sqrt(r1*r1 + r2*r2 -2*r1*r2*cos(theta2)) );
}


double radial_probability_density(double r1,double r2){
  /*
  Normalized radial probability distribution.
  */
  return 16*exp(-4*(r1+r2));
}

double test_func(double x, double y, double z){
  return x*x*y*y*z*z*exp(-(x+y+z));
}



//  Note that you need to call it with a given value of alpha,
// called alf here. This comes from x^{alpha} exp(-x)

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


double laguerre_integrate_func_6(double u1, double u2, double theta1, double theta2, double phi1, double phi2){
  if (u1*u1 + u2*u2 - 2*u1*u2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1 - phi2)) <  1e-8){
    return 0;
  }

  return exp(-3*(u1+u2))*(sin(theta2)*sin(theta1)*u1*u1*u2*u2/sqrt(u1*u1 + u2*u2 - 2*u1*u2*(cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(phi1 - phi2))));

}

double laguerre_integrate_func_3(double u1, double u2, double theta2){
  if(sin(theta2)/sqrt(u1*u1 + u2*u2 - 2*u1*u2*cos(theta2)) < 1e-8){
    return 0;
  }

  return (sin(theta2)/sqrt(u1*u1 + u2*u2 - 2*u1*u2*cos(theta2)));

}
