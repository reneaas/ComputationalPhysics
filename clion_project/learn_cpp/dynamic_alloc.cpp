#include <iostream>
using namespace std;

int main(int argc, char * argv[]){
  int i, npts;
  double *x;      //Declaration of pointer variable 'x'
  double dx;

  cout << "Enter the number of points in [0,1]: ";
  cin >> npts;

  x = new double[npts];         //Dynamic allocation of npts doubles
  dx = 1.0/(npts-1);

  for (int i = 0; i < npts; i++){
    x[i] = i*dx;
  }

  for (int i = 0; i < npts; i++){
    cout << "x [" << i << "] = " << x[i] << endl;
  }

  delete[] x;     //deallocation of dynamically allocated memory.
  return 0;
}
