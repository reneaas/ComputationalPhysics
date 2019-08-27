#include <iostream>

using namespace std;


int main(){
  double a = 1.0;
  double machine_zero = 1.0;
  while (a + machine_zero != 1){
    machine_zero /= 2.0;
  }
  cout << machine_zero << endl;
  return 0;
}
