CXXFLAGS= -std=c++11 -g -Wall  -fopenmp

OBJS = sudoku8tile.o Sudokoid_8tile.o

8tile: sudoku8tile.o Sudokoid_8tile.o
	g++ ${OBJS} $(CXXFLAGS) -o 8tile -lm

sudoku8tile.o: sudoku8tile.cpp Sudokoid_8tile.o 
	g++  $(CXXFLAGS)  -c sudoku8tile.cpp -lm

Sudokoid_8tile.o: Sudokoid_8tile.cpp Sudokoid_8tile.h
	g++  $(CXXFLAGS)  -c Sudokoid_8tile.cpp -lm

clean:
	rm -rf *o 8tile
