from os import system
import sys

number_of_gridpoints = [int(10**(i)) for i in range(1,N+1)]

print("Compiling program...")
system("c++ -O3 -c project1_LU.cpp functions.cpp -larmadillo")
system("c++ -O3 -o project1_LU.exe project1_LU.o functions.o -larmadillo")


print("Compilation done, executing")
for n in number_of_gridpoints:
    print("computing for n = " + str(n))
    system("./project1_LU.exe" + " " + str(n))
