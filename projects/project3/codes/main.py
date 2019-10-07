import os


print("compiling")
os.system("c++ -O3 -c main.cpp lib.cpp")
os.system("c++ -O3 -o main.exe main.cpp lib.o")
print("executing")
os.system("./main.exe")
