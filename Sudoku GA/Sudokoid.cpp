
#include <iostream>
#include <vector>
#include <random>
#include <time.h>
#include <fstream>
#include <exception>
#include <cstdlib>
#include <climits>
#include "SimpleSolver.h"

using namespace std;

/******************************************************************************
copy
copies a given 2-d array of intos into the current cell

Parameters:
	cell - the 2d array to copy from
	cellDim - the number of columns and rows in the 2d array
******************************************************************************/
void Sudokoid::SudokuCell::copy( char **cell, int cellDim)
{
	CellDimension = cellDim;
	Cell.resize(CellDimension);

	for(int i=0; i < CellDimension ; i++){
		Cell[i].resize(CellDimension);
		 for(int j=0; j < CellDimension ; j++){
			Cell[i][j] = cell[i][j];
		  }
		}
}

/******************************************************************************
copy
copies a given vector of 2d ints to copy into the current cell

Parameter:
	cell - the values to copy into the current cell
******************************************************************************/
void Sudokoid::SudokuCell::copy( vector< vector< char > > cell)
{
		Cell = cell;
		CellDimension = cell.size();
}

/******************************************************************************
copy
copies a given sudokucell into the current cell

Parameter:
	sudokuCell - the sudokuCell to copy
******************************************************************************/
void Sudokoid::SudokuCell::copy( SudokuCell sudokuCell)
{
	copy(sudokuCell.Cell);
}

/******************************************************************************
SudokuCell

constructs a new sudoku cell with cell dimension of cellDimension
******************************************************************************/
Sudokoid::SudokuCell::SudokuCell( int cellDimension )
{
	CellDimension = cellDimension;
	Cell.resize(cellDimension);
	for(int i = 0; i < cellDimension; i++)
	{
		Cell[i].resize(cellDimension);
	}
}

/******************************************************************************
SudokuCell

constructs a new sudoku cell assuming a dimension of 3.
Not the ideal solution to some problems, but works for
this implementation.
******************************************************************************/
Sudokoid::SudokuCell::SudokuCell(  )
{
	CellDimension = 3;
	Cell.resize(CellDimension);
	for(int i = 0; i < CellDimension; i++)
	{
		Cell[i].resize(CellDimension);
	}
}

/******************************************************************************
fillBlanks
goes through every character in the cell, replacing '-' characters with
random, unused numbers
******************************************************************************/
void Sudokoid::SudokuCell::fillBlanks( )
{
	int numUsed = 0;
	bool used[9]; //used[n] tells if n already exists in the cell
	fill_n(used, 9, false);
	
	//look through every cell, marking down used numbers
	for(int i=0; i < CellDimension ; i++)
	{
		 for(int j=0; j < CellDimension ; j++)
		 {
			if(Cell[i][j] != '-')
			{
				if (!used[(int)(Cell[i][j] - '0' - 1)])
				{
					used[(int)(Cell[i][j] - '0' - 1)] = true; //convert the char to an int and mark it as used
					numUsed ++;
				}
			}
		}
	}

	//replace blanks with random chars
	for(int i=0; i < CellDimension ; i++)
	{
		for(int j=0; j < CellDimension ; j++)
		{
			if(numUsed > 8)
				return; //there are no blanks to fill
			if (Cell[i][j] == '-')
			{			
				int randnum = -1;
				while(randnum < 0 || used[randnum - 1])
				{
					randnum = (rand() % 9) + 1;
				}
				Cell[i][j] = (char)(((int)'0') + randnum); //convert the random int to a char
				used[randnum - 1] = true; //mark that number as used to prevent repeats
				numUsed++;
			}
			
		}
	}
}

/******************************************************************************
updateLock
goes through every character in the cell, setting the equivalent lock bool
to true if the character is a '-'
******************************************************************************/
void Sudokoid::SudokuCell::updateLock( )
{
	Lock.resize(CellDimension);
	for(int i=0; i < CellDimension ; i++)
	{
		Lock[i].resize(CellDimension);
		 for(int j=0; j < CellDimension ; j++)
		{
			if(Cell[i][j] != '-')
				Lock[i][j] = true;
			else
				Lock[i][j] = false;
		}
	}
}

/******************************************************************************
mutate
swaps two random positions in the sudoku cell
******************************************************************************/
Sudokoid::SudokuCell Sudokoid::SudokuCell::mutate()
{
	//to prevent infinite loops
	if(Locked)
		return *this;	

	int unlockedCount = 0;
	//check to see if there are only 1 or 0 unlocked spots
	for(int i  = 0; i < CellDimension; i++)
		for(int j = 0; j < CellDimension; j++)
			if(!Lock[i][j])
				unlockedCount++;

	//lock the cell to prevent future mutation calls on it.
	if(unlockedCount < 2)
	{
		Locked = true;
		return *this;
	}

	int positionA = rand()%(CellDimension*CellDimension); // get a random position

	//while the first position is locked or not generated yet, set it to a new random position	
	while ( Lock[positionA/CellDimension][positionA%CellDimension])
		positionA = rand()%(CellDimension*CellDimension); 
	
	int positionB = rand()%(CellDimension*CellDimension); //get another random position which is not the same as A

	//while the second position is the same as the first, locked, or not generated yet, set it to a new random position
	while ( Lock[positionB/CellDimension][positionB%CellDimension] || positionB == positionA) 
		positionB = rand()%(CellDimension*CellDimension);

	//swap the two 
	swap(Cell[positionA/CellDimension][positionA%CellDimension], Cell[positionB/CellDimension][positionB%CellDimension]);

	return *this;
}

/******************************************************************************
Sudokoid
a basic constructor for a sudokoid
******************************************************************************/
Sudokoid::Sudokoid()
{
	Fitness = INT_MAX;
	Dimension = 0;
}

/******************************************************************************
Sudokoid
constructs a sudokoid, containing one solution

Parameters:
	values - a pointer to the 2d array of sudoku cells to be used in 
		the solution
	dim - the number of CELLS in the grid
******************************************************************************/
Sudokoid::Sudokoid(SudokuCell **cells, int dimension)
{
	Fitness = INT_MAX;
	Dimension = dimension;
	Cells.resize(Dimension);
		for(int i=0; i < Dimension ; i++){
			Cells[i].resize(Dimension);
			 for(int j=0; j < Dimension ; j++){
				Cells[i][j] = SudokuCell(Dimension);
				Cells[i][j].copy(cells[i][j]);
			  }
			}
}

/******************************************************************************
sudokoid
constructs a sudokoid, containing one solution

Parameter:
	values - a pointer to the 2d array of sudoku cells to be used in the 
	solution
******************************************************************************/
Sudokoid::Sudokoid(vector< vector< SudokuCell > > cells)
{
	Fitness = INT_MAX;
	Dimension = cells.size();
	Cells.resize(Dimension);
		for(int i=0; i < Dimension ; i++){
			Cells[i].resize(Dimension);
			 for(int j=0; j < Dimension ; j++){
				Cells[i][j].copy(cells[i][j]);
			  }
			}
}

/*****************************************************************************
operator<

 Overload of the < operator, allowing use of std::sort
 for sorting vectors of Sudokoids.
 Note - unfitted sudokus will always have INT_MAX fitness.

*****************************************************************************/
bool Sudokoid::operator< (const Sudokoid &other) const 
{
	//compare the two and return the result
	return Fitness < other.Fitness;
}

/******************************************************************************
Fit
Stores the fitness of the current Sudokoid by evaluating how close 
it is to a solution.

100 means a perfect solution, 
while 0 means that no cells were complete and no rows were complete.
******************************************************************************/
void Sudokoid::Fit()
{
	Fitness = SimpleSolver::ssolve(((*puzzle).puzzle.solution));
}

/******************************************************************************
Fit
Stores the fitness of the current Sudokoid by evaluating how close
it is to a solution.

100 means a perfect solution,
while 0 means that no cells were complete and no rows were complete.
******************************************************************************/
int Sudokoid::Fit(Puz &sol)
{
	int fitness = SimpleSolver::ssolve(sol.solution);
	return fitness;
}

/******************************************************************************
Sudokoid
constructs a sudokoid, containing one solution

Parameter:
	puzzle - a puzzle from which to construct the Sudokoid
The number of cells in the grid is assumed to be 9 x 9, since a puzzle is 
	being used
******************************************************************************/
Sudokoid::Sudokoid(Puzzle *puzzleIn)
{
	Dimension = 3;
	puzzle = puzzleIn;
	Fitness = INT_MAX;

	//set up the dimensions
	Cells.resize(Dimension);
	for( int i = 0; i < 3; i++ ) //row
	{
		Cells[i].resize(Dimension);
		for( int j = 0; j < 3; j++ ) //col
		{
			Cells[i][j] = SudokuCell(Dimension);
		}
	}

	for( int i = 0; i < 9; i++ ) //row
	{
		for( int j = 0; j < 9; j++ ) //col
		{
			(Cells[i/3][j/3]).Cell[i%3][j%3] = puzzle->puzzle.solution[i][j][0];
		}
	}
}

/******************************************************************************
Mate
Returns a new offspring Sudokoid which has been crossed over and
mutated to create a new Sudokoid solution.  Crossover is done
on a cellular level, where each cell has the chance of being from
one parent or the other with equal probability.

Mutation is done is done on a subcellular level.  There is a mutationRate
chance any given cell will have a single random swap of two non-locked values.

Parameters:
	mate - the Sudokoid to mate with
	mutationRate - the chance that any given cell is a mutant
******************************************************************************/
Sudokoid Sudokoid::Mate( Sudokoid mate, double mutationRate )
{
	//generate offspring by flipping a coin on every cell
	Sudokoid Sudokling(*this); //start by creating a clone of this
	bool crossover;
	bool mutation;

	//crossover with the mate, then check each cell to see if a mutation is in order
	for(int i=0; i < Dimension ; i++){
			for(int j=0; j < Dimension ; j++){

				//get a random bit: a coin flip
				crossover = (rand() % 2 == 1);
				if (crossover)
					Sudokling.Cells[i][j].copy(mate.Cells[i][j]);

				//get a random number and find if it fits under the integer representation
				//	of the mutation rate
				mutation = bool(rand()%1000000 <= int(mutationRate * 10000));
				if (mutation && !Sudokling.Cells[i][j].Locked)
				{
					Sudokling.Cells[i][j] = Sudokling.Cells[i][j].mutate();
				}
			}
		}
	//create and fit a new puzzle reflecting the child
	Sudokling.puzzle = new Puzzle(&Sudokling);
	DeleteList.push_back(Sudokling.puzzle);
	Sudokling.Fit();
	return Sudokling;
}

/******************************************************************************
lockCells

updates the lock on every cell
******************************************************************************/
void Sudokoid::lockCells()
{
	for(int i=0; i < Dimension ; i++){
		for(int j=0; j < Dimension ; j++){
			Cells[i][j].updateLock();
		}
	}
}

/******************************************************************************
fillCells

Fills blanks in every cell and create a new puzzle with the new random chars
******************************************************************************/
void Sudokoid::fillCells()
{
	for(int i=0; i < Dimension ; i++)
	{
		for(int j=0; j < Dimension ; j++)
		{
			Cells[i][j].fillBlanks();
		}
	}
	puzzle = new Puzzle(this);
	DeleteList.push_back(puzzle);
}

/******************************************************************************
Print

displays the solution in a console output

Parameter:
	cout - a reference to the ostream to pass to puzzle's output function
******************************************************************************/
void Sudokoid::Print( ostream& cout)
{
	puzzle->output(cout);
}

/******************************************************************************
deletePuzzles

delete every puzzle pointer in DeleteList and remove all repeated occurances
of the pointer by traversing backwards.  When a value is erased from a vector,
that index is destroyed, reducing the size of the vector.  Traversing backwards
prevents the skipping of indexes during the deletion.
******************************************************************************/
void Sudokoid::deletePuzzles()
{
	for( int j = 0; j < (int)DeleteList.size(); j++)
	{
		Puzzle* deleted = DeleteList[j];
		delete DeleteList[j];

		//erase all occurances of the address from the vector
		for( int i = DeleteList.size() - 1; i >= 0; i-- ) 
		{
			if( DeleteList[i]==deleted )
				DeleteList.erase( DeleteList.begin() + i );
		}
	}
}
