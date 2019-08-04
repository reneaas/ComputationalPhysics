#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;

double f(double);


int main(int argc, char* argv[]){
	double I = 0;
	double x = 0;
	double dx;
	int N;
	double a, b;								//Integration interval [a,b]
	a = atof(argv[1]);
	b = atof(argv[2]);
	N = atoi(argv[3]);					//Defines the number of intergration points.
	dx = (b-a)/N;								//stepsize.

	//Here comes the actual integration part.
	I += 0.5*(f(a) + f(b));
	for (int i = 1; i < N-1; i++){
		I += f(x + i*dx);
	}
	I = I*dx;
	cout << "Integral = " << I << endl;
	return 0;
}


//Defines the function we want to approximate the integral of.
double f(double x){
	double y = sin(x)/x;
	return y;
}
