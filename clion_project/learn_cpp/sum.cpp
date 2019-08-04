#include <iostream>
using namespace std;

int main(){
  double* s;        //Declares pointer
  int n;
  cout << "Write how many terms in the sum: " << endl;
  cin >> n;
  s = new double[n];
  for (int i = 0; i < n; i++){
    s[i] = 1.0/(( (double) i+1)*((double) i+1));
  }
  double total_sum = 0;

  for (int i = 0; i < n; i++){
    total_sum += s[i];
  }

  delete[] s;
  cout << " = " << total_sum << endl;
  return 0;
}
