CXXFLAGS= -std=c++11 -Wall -O3


OBJS = sudoku.o SimpleSolver.o puzzle.o Sudokoid.o


sudoku: sudoku.o puzzle.o SimpleSolver.o Sudokoid.o

	g++ ${OBJS} ${CXXFLAGS} -lm -o sudoku


sudoku.o: sudoku.cpp SimpleSolver.h 

	g++ ${CXXFLAGS} -c sudoku.cpp


SimpleSolver.o: SimpleSolver.cpp SimpleSolver.h puzzle.h

	g++ ${CXXFLAGS} -c SimpleSolver.cpp


puzzle.o: puzzle.cpp puzzle.h

	g++ ${CXXFLAGS} -c puzzle.cpp


Sudokoid.o: Sudokoid.cpp Sudokoid.h

	g++ ${CXXFLAGS} -c Sudokoid.cpp


clean:

	rm -rf *o sudoku
