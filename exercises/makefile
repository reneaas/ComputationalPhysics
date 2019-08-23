CC = c++ -Wall
PROG= 3_1
LIB_FLAGS = -lblas -llapack -larmadillo 
${PROG} :	${PROG}.o
					${CC} ${PROG}.o -o ${PROG}

${PROG}.o : ${PROG}.cpp
						${CC} -c ${PROG}.cpp