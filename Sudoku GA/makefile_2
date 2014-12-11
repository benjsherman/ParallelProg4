CXXFLAGS= -std=c++11 -Wall -O3
MPICXX = mpic++

OBJS = sudoku.o SimpleSolver.o puzzle.o Sudokoid.o


sudoku: sudoku.o puzzle.o SimpleSolver.o Sudokoid.o

	$(MPICXX) $(OBJS) $(CXXFLAGS) -lm -o sudoku


sudoku.o: sudoku.cpp SimpleSolver.h 

	$(MPICXX) $(CXXFLAGS) -c sudoku.cpp


SimpleSolver.o: SimpleSolver.cpp SimpleSolver.h puzzle.h

	$(MPICXX) $(CXXFLAGS) -c SimpleSolver.cpp


puzzle.o: puzzle.cpp puzzle.h

	$(MPICXX) $(CXXFLAGS) -c puzzle.cpp


Sudokoid.o: Sudokoid.cpp Sudokoid.h

	$(MPICXX) $(CXXFLAGS) -c Sudokoid.cpp


clean:

	rm -rf *o sudoku
