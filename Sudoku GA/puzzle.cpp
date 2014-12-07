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
			_in >> puzzle.solution[i][j][0];
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
			puzzle.solution[i][j][0] = (sudokoid->Cells[i/3][j/3]).Cell[i%3][j%3];
		}
	}
}

/******************************************************************************
Puzzle Copy Constructor
Creates a puzzle from an existing puzzle
Parameter:
	Puzzle* const puzzle - puzzle to copy values from
******************************************************************************/
Puzzle::Puzzle(Puzzle* const puz)
{
	for( int i = 0; i < 9; i++ ) //row
	{
		for( int j = 0; j < 9; j++ ) //col
		{
			puzzle.solution[i][j][0] = puz->puzzle.solution[i][j][0];
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
			_out << puzzle.solution[i][j][0] << " ";
			if(j == 2||j == 5 || j == 8)
				_out << " ";
			
		}
		_out << endl;
		if(i == 2 || i == 5 || i == 8)
			_out << endl;
	}
}		
