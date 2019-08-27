#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

ofstream ofile;

void differentiate1(double, double, double&);
void differentiate2(double, double, double&);
void write_to_file(double*, double*,double*, char*);
void compute_error(double, double&);
void differentiate3(double, double, double&);
double f(double);


int main(int argc, char* argv[]){
  char * outfilename;
  int n = 50;
  outfilename = argv[1];
  double *h, *derivative, *error;
  h = new double[n];
  error = new double[n];
  derivative = new double[n];
  double x = M_PI_2;

  //Fill h with a range of stepsizes.
  for (int i = 0; i < n; i++){
    h[i] = pow(0.5, (double) i);
  }

  //Here we compute the derivative.
  for (int i = 0; i < n; i++){
    differentiate3(x, h[i], derivative[i]);
    compute_error(derivative[i], error[i]);
  }


  write_to_file(h, derivative, error, outfilename);
  return 0;
}


double f(double x){
  return sin(x);
}

void compute_error(double derivative, double& error){
  error = log10(abs((-1.0 - derivative)/(-1.0)));
}

void differentiate1(double x, double h, double& derivative){
  derivative = (atan(x+h) - atan(x))/h;
  return;
}

void differentiate2(double x, double h, double& derivative){
  derivative = (atan(x+h) - atan(x-h))/(2*h);
  return;
}

void differentiate3(double x, double h, double& derivative){
  derivative = (f(x+h) - 2*f(x) + f(x-h))/(h*h);
}

void write_to_file(double* h, double* derivative, double* error, char* outfilename){
  ofile.open(outfilename);  //Opens the file we'll write the data to.
  ofile << "h" << " "<< "derivative" << endl;
  int n = 50;
  for (int i = 0; i < n; i++){
    ofile << h[i] << " " << derivative[i] << " " << error[i] << " " << endl;
  }
  ofile.close();
  return;
}
