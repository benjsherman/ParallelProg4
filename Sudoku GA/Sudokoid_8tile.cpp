#include <omp.h>
#include <iostream>
#include <vector>
#include <random>
#include <time.h>
#include <fstream>
#include <exception>
#include <cstdlib>
#include <climits>
#include "Sudokoid_8tile.h"

using namespace std;


//debugging functions

void print_moves(Sudokoid sudokoid)
{  
   cout << std::endl;
   for(int j = 0; j < sudokoid.Moves.Length; j++)
      std::cout << sudokoid.Moves.Moves[j];
   cout << std::endl;
}


/**************************************/
/*	Some non-member utility functions */
/**************************************/


/*****************************************************************************
printGrid

Prints a single 3x3 8tile grid to the screen.

Paramaters:
	grid - the 3x3 grid to print
	cout - the stream to print to

*****************************************************************************/
void printGrid(int grid[][3], ostream& cout)
{
    cout <<"\n\n";
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            cout << grid[i][j] << " ";
        }
        cout << endl;
    }
}

/*****************************************************************************
getScore

 gets the Score of a single 3 x 3 8tile grid and returns it to the 
 fitness function.  Utilizes the global variable Puzzle.
 
 Parameter grid - the grid to evaluate.

*****************************************************************************/
int getScore(int grid[][3])
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




// the starting puzzle, a 3x3
int Puzzle[3][3];

/*****************************************************************************
readPuzzle

Reads in the global 3x3 Puzzle from a specified stream.
The first two values are the width and height, and the rest of the values are
the tiles numbers at their positions.

Paramater istream - the stream to read from.

*****************************************************************************/
void readPuzzle(istream & _in)
{
	int dim1, dim2;
	
	_in >> dim1;
	_in >> dim2;
	if(dim1 != 3 || dim2 != 3)
	{
		throw new string("Only 3 by 3 puzzles are currently supported");
	}
	
	for(int i=0; i < dim1 ; i++)
	{
		for(int j=0; j < dim2; j++)
		{
			_in >> Puzzle[i][j];
		}
	}
}

/*****************************************************************************
printPuzzle

Prints a single 3x3 8tile grid to the screen.

Paramaters:
	grid - the 3x3 grid to print
	cout - the stream to print to

*****************************************************************************/
void printPuzzle(ostream& cout)
{
    printGrid(Puzzle, cout);
}


/*****************************/
/* Sudokoid Member Functions */
/*****************************/


/*****************************************************************************
operator<

 Overload of the < operator, allowing use of std::sort
 for sorting vectors of Sudokoids.
 Note - unfitted sudokus will always have INT_MAX fitness.

*****************************************************************************/
bool Sudokoid::operator< (const Sudokoid &other) const 
{
	//compare the two and return the result
	return this->Fitness < other.Fitness;
}

//constructors

Sudokoid::Sudokoid()
{
   this->Fitness = INT_MAX;
   this->Moves.Fitness = INT_MAX;
   this->Moves.Length = 100;
}

Sudokoid::Sudokoid(const Sudokoid& other)
{
	//copy constructor
    this->copySolution(other.Moves);
    this->Fitness = other.Fitness;
}


/*****************************************************************************
generateRandomSolution

Fills all 100 spots in the solution with random moves

*****************************************************************************/
void Sudokoid::generateRandomSolution()
{
    for(int i=0; i < 100; i++)
    {
        switch(rand()%4)
        {
            case 0: this->Moves.Moves[i] = 'u';
                break;
            case 1: this->Moves.Moves[i] = 'd';
                break;
            case 2: this->Moves.Moves[i] = 'l';
                break;
            case 3: this->Moves.Moves[i] = 'r';
                break;
            default:
               break;
        }
    }
    this->Moves.Length = 100;
}

/*****************************************************************************
copySolution

Copies all elements from one solution to another.

Parameter solution - the solution to copy
*****************************************************************************/
void Sudokoid::copySolution( const Solution solution )
{
	for(int i = 0; i < solution.Length; i++)
	{
		this->Moves.Moves[i] = solution.Moves[i];
	}
	this->Moves.Length = solution.Length;
	this->Moves.Fitness = solution.Fitness;
	this->Fitness = this->Moves.Fitness;
}


/*****************************************************************************
Fit

Steps through the puzzle steps and assigns a fitness to the solution

*****************************************************************************/
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


    if(debugging)
    {
       cout << "\nFitting Sudokoid of Length " << (*this).Moves.Length << ": ";
       print_moves(*this);
       cout << endl;
    }

    for (int i = 0; i < (*this).Moves.Length; i++)
    {
        if(debugging)
        {
            cout << "\nMove " << i << ": " << (*this).Moves.Moves[i] << endl;
            //cout << "\n Zero Row: " << zero_row<< endl;
            //cout << "\n Zero Col: " << zero_col<< endl;
            //printGrid(puzzle, std::cout);
            cout << endl;
        }
        bool valid = false;
        switch( (*this).Moves.Moves[i] )
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
            if(zero_col > 0)
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row][zero_col-1];
                puzzle[zero_row][zero_col-1] = 0;
                zero_col -= 1;
                valid = true;
            }
            break;
        default:
            break;
        }

        if(valid)
        {
            int gotten_score = getScore(puzzle);
            if(gotten_score < best_score)
            {
               best_score = gotten_score;
               best_index = i;

               if(debugging)
               {
                  cout << "New Best!  Score of " << best_score << " at index " << best_index << endl;
                  printGrid(puzzle, std::cout);
               } 
            }
        }
    }
    if(debugging)
    {
       cout << "\n----------------------------------------\nFitness of " << best_score * 100 + best_index;
       cout << endl;
    }

    this->Moves.Length = best_index + 1;
    this->Moves.Fitness = best_score * 100 + best_index;
	 this->Fitness = this->Moves.Fitness;
}

/*****************************************************************************
Print

Prints the solution steps and the 8tile puzzle's state at each step.

Parameter cout - the stream to print to

*****************************************************************************/
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
            if(zero_col > 0)
            {            
                puzzle[zero_row][zero_col] = puzzle[zero_row][zero_col-1];
                puzzle[zero_row][zero_col-1] = 0;
                zero_col -= 1;
                valid = true;
            }
            break;
        default:
            break;
        }

        if(valid)
        {
            printGrid(puzzle, cout);
        }
    }
}

/******************************************************************************
Mate

Returns a new offspring Sudokoid which has been crossed over and
mutated to create a new Sudokoid solution.

Parameters:
	mate - the Sudokoid to mate with
	mutationRate - the chance that any given move is a mutant
******************************************************************************/
Sudokoid Sudokoid::Mate( Sudokoid mate, double mutationRate )
{
	//generate offspring by flipping a coin on every cell
	 Sudokoid Sudokling;
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
                        break;
                }
				Sudokling.Moves.Moves[i] = mutator;
			}
	}

   Sudokling.Fit();

	return Sudokling;
}

