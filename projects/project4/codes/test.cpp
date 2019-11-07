#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstdlib>
#include <string>
#include "time.h"
#include <random>
#include <fstream>
using namespace std;


int main(int nargs, char* args[]){
  int n = atoi(args[1]);
  double a = 0;
  for (long long int i = 0; i < n; i++){
    a += 2;
  }
  cout << a << endl;
}
