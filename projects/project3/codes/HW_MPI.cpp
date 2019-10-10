
#include "mpi.h"
#include <iostream>
using namespace std;
int main(int argc, char** argv)
{

int my_PE_num;
MPI::Init(argc, argv);
my_PE_num = MPI::COMM_WORLD.Get_rank();
cout << "Hello from " << my_PE_num << endl;
MPI::Finalize();
}
