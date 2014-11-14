#include"SimpleSolver.h"

using namespace std;

/******************************************************************************
Simple Solver
This function solves sudoku puzzles by finding naked singles and hidden 
singles. It will return when it can no longer find an values by using these
statagies even if the puzzle is not completly solved. It uses the third 
demintion of the Puzzle->solution array to store possiblities.
Parameter:
	Puzzle &p - an instance of the Puzzle class containing the puzzle
		to be solved
******************************************************************************/
void SimpleSolver::ssolve(Puzzle &p)
{
	int score, score2;
	
	do
	{
		score = p.getSolutionScore();
		
		resetPossiblities(p);
		
		mrkRowColmPossiblities(p);
		
		//narrow possiblities in squares
		markSqrPossiblities(0, 0, p);
		markSqrPossiblities(0, 3, p);
		markSqrPossiblities(0, 6, p);
		markSqrPossiblities(3, 0, p);
		markSqrPossiblities(3, 3, p);
		markSqrPossiblities(3, 6, p);
		markSqrPossiblities(6, 0, p);
		markSqrPossiblities(6, 3, p);
		markSqrPossiblities(6, 6, p);
		
		findNakedSingles(p);
		
		findHiddenSingleRowColn(p);
		
		//check squares for hidden singles
		hiddenSingleSquare(0, 0, p);
		hiddenSingleSquare(0, 3, p);
		hiddenSingleSquare(0, 6, p);
		hiddenSingleSquare(3, 0, p);
		hiddenSingleSquare(3, 3, p);
		hiddenSingleSquare(3, 6, p);
		hiddenSingleSquare(6, 0, p);
		hiddenSingleSquare(6, 3, p);
		hiddenSingleSquare(6, 6, p);			 
		
		//return if no progress is being made
		score2 = p.getSolutionScore();	
		if(score2 == score)
			return;
					
	}while(score > 0);
}

/******************************************************************************
Mark Square Possiblities
Marks the impossible values as '0' in a 3 x 3 square starting with column x
and row y.
Parameters:
	Puzzle &p - the puzzle
	int x - starting location in the 1st diminsion
	int y - starting location in the 2nd diminsion
******************************************************************************/ 
void SimpleSolver::markSqrPossiblities(int x, int y, Puzzle &p)
{
	int i, j, m, k = 0;
	char tarry[9];
	
	for(i = x; i < x+3; i++)
	{
		for(j = y; j < y+3; j++)
		{
			tarry[k] = p.solution[i][j][0];
			k++;
		}
	}
	for(m = 0; m < 9; m++)
	{
		if(tarry[m] != '-')
		{
			k = tarry[m] - '0';
			for(i = x; i < x+3; i++)
			{
				for(j = y; j < y+3; j++)
				{
					p.solution[i][j][k] = '0';	
				}
			}
		}
	}
}

/******************************************************************************
Hidden Singles Square
Finds the hidden singles in 3 x 3 squares
Parameters:
	int x - starting location in the 1st deminsion
	int y - starting location in the 2nd deminsion
	Puzzle &p - the puzzle
******************************************************************************/
void SimpleSolver::hiddenSingleSquare(int x, int y, Puzzle &p)
{
	int i, j, k;
	int n = -1;
	int m = -1;
	

	for(k = 1; k < 10; k++)
	{
		int found = 0;
		for(i = x; i < x+3; i++)
		{
			for(j = y; j < y+3; j++)
			{
				if(p.solution[i][j][k] == '1')
				{
						n = i;
						m = j;
						found++;
				}
			}
		}
		//cout << "found: " << k  << " " << found << " times "<< " in square " << x << " " << y << endl;
		if(found == 1)
		{
			p.solution[n][m][0] = k + '0';
		}
	}
}

/******************************************************************************
Reset Possiblities
Sets possiblilies behind filled places to '0' and possiblities behind blanks to
'1'
Parameter:
	Puzzle &p - the puzzle
******************************************************************************/
void SimpleSolver::resetPossiblities(Puzzle &p)
{
	int i, j, k; 
	char temp;
	//initial setup of the solver
	for(i = 0; i < 9; i++)
	{
		for(j = 0; j < 9; j++)
		{
			if (p.solution[i][j][0] == '-')
				temp = '1';
			else
				temp = '0';
			for(k = 1; k < 10; k++)
				p.solution[i][j][k] = temp;
		}	
	}
}

/******************************************************************************
Mark Row & Column Possiblities
Narrows down the possiblities by checking the rows and columns for values that
would create duplicates
Parameter:
	Puzzle &p - the puzzle
******************************************************************************/
void SimpleSolver::mrkRowColmPossiblities(Puzzle &p)
{
	int i, j, k, x, y;
	//narrow possiblities in rows and columns
	for (i = 0; i < 9; i++)
	{
		for(j = 0; j < 9 ; j++)
		{
			if(p.solution[i][j][0] != '-')
			{
				k = p.solution[i][j][0] - '0';
				for(x = 0; x < 9; x++)
				{
					p.solution[x][j][k] = '0';
				}
				for(y = 0; y < 9; y++)
				{
					p.solution[i][y][k] = '0';
				}
			}
		}
	}
}

/******************************************************************************
Find Naked Singles
Fills in places with only one possiblity marked behind them.
Parameter:
	Puzzle &p - the puzzle
******************************************************************************/
void SimpleSolver::findNakedSingles(Puzzle &p)
{
	int i, j, k, found;;
	//fill in naked singles
	for (i = 0; i < 9; i++)
	{
		for(j = 0; j < 9 ; j++)
		{
			if(p.solution[i][j][0] == '-')
			{
				found = 0;
				for (k = 1; k < 10; k++)
				{
					if(p.solution[i][j][k] == '1')
					{
						if(found == 0)
						{
							found = k;
						}
						else
						{	
							found = -1;
						}
					}
				}

				if(found > 0)
					p.solution[i][j][0] = found + '0';
			}
		}
	}
}

/******************************************************************************
Find Hidden Singles In Rows and Columns
Fills in values which are the only possible place to put that digit in their
row or column.
Parameter:
	Puzzle &p - the puzzle
******************************************************************************/
void SimpleSolver::findHiddenSingleRowColn(Puzzle &p)
{
	int i, j, k, x, y;
	int found;
	//find hidden singles
	for(k = 1; k < 10; k++)
	{
		//check columns
		for(i = 0; i < 9; i++)
		{
			x = -1;
			y = -1;
			found = 0;
			for(j = 0; j < 9; j++)
			{
				if(p.solution[i][j][k] == '1')
				{
					x = i;
					y = j;
					found++;
				}
			}
			if(found == 1)
			{
				p.solution[x][y][0] = k + '0';					
			}
		}
		//check rows
		for(j = 0; j < 9; j++)
		{
			x = -1;
			y = -1;
			found = 0;
			for(i = 0; i < 9; i++)
			{
				if(p.solution[i][j][k] == '1')
				{
					x = i;
					y = j;
					found++;
				}
			}
			if(found == 1)
			{
				p.solution[x][y][0] = k + '0';
			}
		}
	}	
}
