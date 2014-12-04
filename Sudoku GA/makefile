CXXFLAGS= -std=c++11 -g -Wall -O3 -fopenmp


OBJS = sudoku.o SimpleSolver.o puzzle.o Sudokoid.o


sudoku: sudoku.o puzzle.o SimpleSolver.o Sudokoid.o

	g++ ${OBJS} ${CXXFLAGS} -o sudoku -lm


sudoku.o: sudoku.cpp SimpleSolver.h 

	g++ ${CXXFLAGS} -c sudoku.cpp -lm


SimpleSolver.o: SimpleSolver.cpp SimpleSolver.h puzzle.h

	g++ ${CXXFLAGS} -c SimpleSolver.cpp -lm


puzzle.o: puzzle.cpp puzzle.h

	g++ ${CXXFLAGS} -c puzzle.cpp -lm


Sudokoid.o: Sudokoid.cpp Sudokoid.h

	g++ ${CXXFLAGS} -c Sudokoid.cpp -lm


clean:

	rm -rf *o sudoku
