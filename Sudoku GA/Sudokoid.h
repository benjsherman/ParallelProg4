#include <omp.h>
#include <iostream>
#include <vector>
#include <random>
#include <time.h>
#include <fstream>
#include <exception>
#include <cstdlib>

using namespace std;

class Puzzle; //forward declare the Puzzle file, which we couldn't do because of circular includes

/******************************************************************************
Sudokoid
contains information on a single sudoku solution
for use in genetic algorithm functions in sudoku.cpp

author - Derek Stotz

SudokuCell
holds information in a single cell of the sudoku grid 
(a 9x9 grid would have 9 of these)
******************************************************************************/
class Sudokoid
{

	public:
		struct SudokuCell
		{
			int CellDimension;	//The width/height of the sudoku cell
			vector< vector< char > > Cell;	//The values in the cell, organized in a 2-d array
			vector< vector< bool > > Lock;	//Whether or not each cell is locked, per the solution
			bool Locked; //Whether or not the entire SudokuCell is locked to mutation (only 1 or 0 unlocked numbers)
	
			//construcutors
			SudokuCell( int cellDimension);
			SudokuCell();

			void copy( char **cell, int cellDim);
			void copy( vector< vector< char > > cell);
			void copy( SudokuCell sudokuCell);
			SudokuCell mutate();
			void fillBlanks();
			void updateLock();
		};

		vector < vector < SudokuCell > > Cells; //the Dimension by Dimension cells in the grid
		int Dimension; //the number of cells wide the grid is
		int Fitness; //the fitness of the current sudokoid
		Puzzle *puzzle; //should be a Puzzle pointer

		static vector <Puzzle*> DeleteList; //a list of all puzzles on the heap

		//constructors
		Sudokoid();
		Sudokoid(SudokuCell **cells, int dimension);
		Sudokoid(Puzzle *puzzle);
		Sudokoid(vector< vector< SudokuCell > > cells);
		
		//other methods
		bool operator< (const Sudokoid &other) const;
		void fillCells();
		void lockCells();
		void Fit();
		Sudokoid Mate( Sudokoid mate, double mutationRate );
		void Print( ostream& cout);
		
		//memory management
		static void deletePuzzles();
};
