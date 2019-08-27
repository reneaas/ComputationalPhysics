#include <iostream>

using namespace std;

void lol(double*, int);

int main(){
  int n = 10;
  double* a = new double[n];
  lol(a,n);
  for (int i =0; i < n; i++){
    cout << l[i] << endl;
  }
  return 0;
}


void lol(double* a, int n){
  double *l = new double[n];
  
  for (int i = 0; i < n; i++){
    a[i] = 1.0;
    l[i] = 1.0;
  }
  return;
}
