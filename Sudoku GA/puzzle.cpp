#include "puzzle.h"
#include <string>

using namespace std;

/******************************************************************************
Puzzle Constructor
Reads in the puzzle from a input stream
Parameter:
	istream &_in - input stream to read puzzle from
******************************************************************************/ 	
Puzzle::Puzzle(istream & _in)
{
	int dim1, dim2;
	
	_in >> dim1 >> dim2;
	if(dim1 != 9 || dim2 != 9)
	{
		throw new string("Only nine by nine puzzles are currently supported");
	}
	
	for(int i=0; i < 9 ; i++)
	{
		for(int j=0; j < 9; j++)
		{
			_in >> solution.sol[i][j][0];
		}
	}
	
}

/******************************************************************************
Sudokoid to Puzzle Constructor
Creates a puzzle from a sudokoid
Parameter:
	Sudokoid* const sudokoid - sudokoid to convert to puzzle
******************************************************************************/
Puzzle::Puzzle(Sudokoid* const sudokoid)
{
	for( int i = 0; i < 9; i++ ) //row
	{
		for( int j = 0; j < 9; j++ ) //col
		{
			solution.sol[i][j][0] = (sudokoid->Cells[i/3][j/3]).Cell[i%3][j%3];
		}
	}
}

/******************************************************************************
Puzzle Copy Constructor
Creates a puzzle from an existing puzzle
Parameter:
	Puzzle* const puzzle - puzzle to copy values from
******************************************************************************/
Puzzle::Puzzle(Puzzle* const puzzle)
{
	for( int i = 0; i < 9; i++ ) //row
	{
		for( int j = 0; j < 9; j++ ) //col
		{
			solution.sol[i][j][0] = puzzle->solution.sol[i][j][0];
		}
	}
}

/******************************************************************************
Output
Prints the puzzle to the give output stream
Paramenter:
	ostream &_out - stream to print to
******************************************************************************/
void Puzzle::output(ostream & _out)
{
	for(int i=0; i < 9 ; i++)
	{
		for(int j=0; j < 9; j++)
		{
			_out << solution.sol[i][j][0] << " ";
			if(j == 2||j == 5 || j == 8)
				_out << " ";
			
		}
		_out << endl;
		if(i == 2 || i == 5 || i == 8)
			_out << endl;
	}
}		

/******************************************************************************
Get solution.sol Score
Calculates a score for the solution.sol taking into account blank spaces and 
duplicates in rows, columns, and squares. Better solution.sols will get lower 
scores. A perfect solution.sol recieves a score of 0.
******************************************************************************/ 
int Puzzle::getSolutionScore()
{
	int score = 0;
	
	//check for blanks
	score = countBlanks();
	
	//check rows for duplicates
	score += rowDuplicates();
	
	//check columns for duplicates
	score += colDuplicates();
	
	//check for duplicates in squares
	score += squareDuplicates(0,0);
	score += squareDuplicates(0,3);	
	score += squareDuplicates(0,6);
	score += squareDuplicates(3,0);
	score += squareDuplicates(3,3);
	score += squareDuplicates(3,6);
	score += squareDuplicates(6,0);
	score += squareDuplicates(6,3);
	score += squareDuplicates(6,6);
 	
	return score;
}	

/******************************************************************************
Count Blanks
Counts the blank spaces in a solution.sol
******************************************************************************/
int Puzzle::countBlanks()
{
	int i, j, score = 0;
	for(i=0; i<9; i++)
	{
		for(j=0; j<9;j++)
		{
			if(solution.sol[i][j][0] == '-')
				score++;
		}
	}
	return score;
}

/******************************************************************************
Row Duplicates
Counts pairs of row duplicates
******************************************************************************/
int Puzzle::rowDuplicates()
{	
	int i, j, k, score = 0;
	char temp;
	for(i = 0; i < 9; i++)
	{
		for(j = 0; j < 9;j++)
		{
			if (solution.sol[i][j][0] != '-')
			{
				temp = solution.sol[i][j][0];
				for(k = j+1; k < 9; k++)
				{
					if (solution.sol[i][k][0] == temp)
					{
						score++;
					}
				}
			}
		}		
	}
	return score;
}

/******************************************************************************
Column Duplicates
Counts Pairs of column duplicates
******************************************************************************/
int Puzzle::colDuplicates()
{
	int i, j, k, score = 0;
	char temp;
	for(j = 0; j < 9 ; j++)
	{
		for(i = 0; i < 9; i++)
		{
			if (solution.sol[i][j][0] != '-')
			{
				temp = solution.sol[i][j][0];
				for(k = i+1;k < 9; k++)
				{
					if (solution.sol[k][j][0] == temp)
					{
						score++;
					}
				}
			}
		}		
	}
	return score;
}

/******************************************************************************
Square  Duplicates
Counts pairs of duplicates in a 3 x 3 square.
Parameters:
	int x - the start of the square in the 1st dimension
	int y - the start of the square in the 2nd dimension
******************************************************************************/
int Puzzle::squareDuplicates(int x, int y)
{
	int i, j, k = 0, score = 0;
	char tarry[9];
	
	for(i = x; i < x+3; i++)
	{
		for(j = y; j < y+3; j++)
		{
			tarry[k] = solution.sol[i][j][0];
			k++;
		}
	}
	for(i = 0; i < 9; i++)
	{
		if(tarry[i] != '-')
		{
			for(k = i+1; k < 9; k++)
			{
				if(tarry[i] == tarry[k])
					score++;
			}
		}
	}
	return score;
}
