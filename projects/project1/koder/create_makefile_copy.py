#Automatically creates a makefile for a c++ program.
print("Type name of your c++ file (without extension):")
PROG = str(input())

with open("makefile", "w") as outfile:
    outfile.write("CC = c++ -Wall" + "\n")
    outfile.write("PROG= " + PROG + "\n")
    outfile.write("LIB_FLAGS = -lblas -llapack -larmadillo ")
    outfile.write("\n")
    outfile.write("${PROG} :	${PROG}.o" + "\n")
    outfile.write("					${CC} ${PROG}.o -o ${PROG}" + "\n")
    outfile.write("\n")
    outfile.write("${PROG}.o : ${PROG}.cpp" + "\n")
    outfile.write("						${CC} -c ${PROG}.cpp")
