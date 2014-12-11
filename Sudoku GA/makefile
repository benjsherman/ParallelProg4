OBJS = sudoku8tile.o Sudokoid_8tile.o

8tile: sudoku8tile.o Sudokoid_8tile.o
	g++ -g ${OBJS} -lm -std=c++11 -o 8tile

sudoku8tile.o: sudoku8tile.cpp Sudokoid_8tile.o 
	g++ -Wall -g -std=c++11 -c sudoku8tile.cpp

Sudokoid_8tile.o: Sudokoid_8tile.cpp Sudokoid_8tile.h
	g++ -Wall -g -std=c++11 -c Sudokoid_8tile.cpp

clean:
	rm -rf *o 8tile
