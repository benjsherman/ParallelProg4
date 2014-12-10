
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

// the starting puzzle, a 3x3
int Puzzle[3][3];

void readPuzzle(istream & _in)
{
	int dim1, dim2;
	
	_in >> dim1 >> dim2;
	if(dim1 != 3 || dim2 != 3)
	{
		throw new string("Only 3 by 3 puzzles are currently supported");
	}
	
	for(int i=0; i < 3 ; i++)
	{
		for(int j=0; j < 3; j++)
		{
			_in >> Puzzle[i][j];
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
	return this->Moves.Fitness < other.Moves.Fitness;
}

//constructors

Sudokoid::Sudokoid()
{

}

Sudokoid::Sudokoid(Solution solution)
{
    this->Moves = solution;
}

void printGrid(int[][3] grid)
{
    std::cout <<"\n\n";
    for (int i = 0; i < 3; i++)
    {
        std::cout <<"\n";
        for (int j = 0; j < 3; j++)
        {
            std::cout << grid[i][j] << " ";
        }
    }
}

//Solution methods
void Sudokoid::copySolution( const Solution solution );

void Sudokoid::generateRandomSolution()
{
    for(int i=0; i < 100; i++)
    {
        switch(rand()%4)
        {
            case 0: this->Moves->Moves[i] = 'u';
                break;
            case 1: this->Moves->Moves[i] = 'd';
                break;
            case 2: this->Moves->Moves[i] = 'l';
                break;
            case 3: this->Moves->Moves[i] = 'r';
                break;
            default:
        }
    }
    this->Moves->Length = 100;
}

void Sudokoid::Fit()
{
    int best_score = 10;
    int best_index = 0;

    // to keep track of the blank
    int zero_row;
    int zero_col;

    // copy the puzzle into a local array
    int puzzle[3][3];
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(Puzzle[i][j] == 0)
            {
                zero_row = i;
                zero_col = j;
            }
            puzzle[i][j] = Puzzle[i][j];
        }
    }

    // Solution looks like
    // 1 2 3
    // 4 5 6
    // 7 8 0

    // start moving

    for (int i = 0; i < this->Moves.Length; i++)
    {
        bool valid = false;
        switch( this->Moves.Moves[i] )
        {
        case 'u':
            if(zero_row > 0)
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row - 1][zero_col];
                puzzle[zero_row-1][zero_col] = 0;
                zero_row -= 1;
                valid = true;
            }
            break;
        case 'd':
            if(zero_row < 2)	
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row + 1][zero_col];
                puzzle[zero_row+1][zero_col] = 0;
                zero_row += 1;
                valid = true;
            }
            break;
        case 'r':
            if(zero_col < 2)
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row][zero_col+1];
                puzzle[zero_row][zero_col+1] = 0;
                zero_col += 1;
                valid = true;
            }
            break;
        case 'l':
            if(zero_row > 0)
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row][zero_col-1];
                puzzle[zero_row][zero_col-1] = 0;
                zero_col -= 1;
                valid = true;
            }
            break;
        default:
        }

        if(valid)
        {
            new_solution.Moves.Moves[new_solution.Length] = this->Moves.Moves[i];
            new_solution.Length++;

            int gotten_score = getScore(puzzle);
            if(gotten_score < best_score)
            {
                best_score = gotten_score;
                best_index = i;
            }
        }
    }

    this->Moves.Length = best //there are no blanks to fill
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
		}_index+1;
    this->Moves.Fitness = best_score * 100 + best_index;
}

int getScore(int[][3] grid)
{
    int score = 0;

    for (int i = 1; i < 9; i++)
    {
        if(grid[(i-1)/3][(i-1)%3] != i)
            score++;
    }

    if (grid[2][2] != 0)
        score++;

    return score;
}

void Sudokoid::Print( ostream& cout)
{
    // to keep track of the blank
    int zero_row;
    int zero_col;

    // copy the puzzle into a local array
    int puzzle[3][3];
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            if(Puzzle[i][j] == 0)
            {
                zero_row = i;
                zero_col = j;
            }
            puzzle[i][j] = Puzzle[i][j];
        }
    }

    // Solution looks like
    // 1 2 3
    // 4 5 6
    // 7 8 0

    // start moving

    for (int i = 0; i < this->Moves.Length; i++)
    {
        bool valid = false;   
        switch( this->Moves.Moves[i] )
        {
        case 'u':
            if(zero_row > 0)
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row - 1][zero_col];
                puzzle[zero_row-1][zero_col] = 0;
                zero_row -= 1;
                valid = true;
            }
            break;
        case 'd':
            if(zero_row < 2)	
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row + 1][zero_col];
                puzzle[zero_row+1][zero_col] = 0;
                zero_row += 1;
                valid = true;
            }
            break;
        case 'r':
            if(zero_col < 2)
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row][zero_col+1];
                puzzle[zero_row][zero_col+1] = 0;
                zero_col += 1;
                valid = true;
            }
            break;
        case 'l':
            if(zero_row > 0)
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row][zero_col-1];
                puzzle[zero_row][zero_col-1] = 0;
                zero_col -= 1;
                valid = true;
            }
            break;
        default:
        }

        if(valid)
        {
            printGrid(puzzle);
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
	Sudokoid Sudokling();
    int cross_index;
    int mate_start, mate_end, this_start, this_end;

    if(this->Moves.Length < mate.Moves.Length)
    {
        cross_index = rand() % this->Moves.Length;
        Sudokling.Moves.Length = this->Moves.Length;

        // I have the faster solution!  I should be the second part.
        mate_start = 0;             mate_end = cross_index;
        this_start = cross_index;   this_end = this->Moves.Length;
    }
    else
    {
        cross_index = rand() % mate.Moves.Length;
        Sudokling.Moves.Length = mate.Moves.Length;

        this_start = 0;             this_end = cross_index;
        mate_start = cross_index;   mate_end = this->Moves.Length;
    }

    for(int i=mate_start; i < mate_end; i++)
        Sudokling.Moves.Moves[i] = mate.Moves.Moves[i];
    for(int i=this_start; i < this_end; i++)
        Sudokling.Moves.Moves[i] = this->Moves.Moves[i];
    
    // apply mutation
    for(int i=0; i<mate.Moves.Length; i++)
    {
			//get a random number and find if it fits under the integer representation
			//	of the mutation rate
			if (rand()%1000000 <= int(mutationRate * 10000))
			{
                char mutator;
                switch(rand()%4)
                {
                    case 0: mutator = 'u';
                        break;
                    case 1: mutator = 'd';
                        break;
                    case 2: mutator = 'l';
                        break;
                    case 3: mutator = 'r';
                        break;
                    default:
                }
				Sudokling.Moves.Moves[i] = mutator;
			}
	}

	Sudokling.Fit();
	return Sudokling;
}

