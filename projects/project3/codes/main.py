import os
import sys

compilation_instruction = str(sys.argv[1])
if compilation_instruction == "no_mpi":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    print("executing")
    os.system("./main.exe")


if compilation_instruction == "0":
    print("executing")
    os.system("./main.exe")

if compilation_instruction == "mpi":
    print("compiling with mpi")
    os.system("mpicxx -O3 -c main_mpi.cpp")
    os.system("mpicxx -O3 -o main_mpi.exe main_mpi.o")
    print("executing with mpi")
    os.system("mpirun -np 10 ./main_mpi.exe")
