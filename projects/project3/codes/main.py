import os
import sys

compilation_instruction = str(sys.argv[1])
if compilation_instruction == "1":
    print("compiling")
    os.system("c++ -O3 -c main.cpp lib.cpp")
    os.system("c++ -O3 -o main.exe main.cpp lib.o")
    print("executing")
    os.system("./main.exe")


if compilation_instruction == "0":
    print("executing")
    os.system("./main.exe")
