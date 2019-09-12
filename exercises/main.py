from os import system
import sys

n = int(sys.argv[1])


system("c++ -O3 -c main.cpp")
system("c++ -O3 -o main.exe main.o")

system("./main.exe" + " " + str(n))
