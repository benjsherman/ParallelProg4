#include <mpi.h>
#include <iostream>
#include <vector>
#include <random>
#include <time.h>
#include <fstream>
#include <exception>
#include <cstdlib>

using namespace std;

/* User defined Macros */
#define BLOCK_LOW(id,p,n)	((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)	(BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n)	(BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n)	(((p)*((j)+1)-1)/(n))

// to read in the global puzzle
void readPuzzle(istream &_in);
void printPuzzle(ostream& cout);
void printGrid(int[][3] grid, ostream& cout);

/******************************************************************************
Sudokoid
contains information on a single 8tile puzzle solution
for use in genetic algorithm functions in sudoku.cpp

author - Derek Stotz

Solution
holds information for the moves made by the blank space in the puzzle
******************************************************************************/
class Sudokoid
{

	public:
		struct Solution
		{
            // cap the solution length at 100

			int Length;	//The number of moves in the solution.
            int Fitness = INT_MAX;
			char Moves[100];	//The moves
		};

		int Fitness = INT_MAX;
		Solution Moves;  //the solution data

		//constructors
		Sudokoid();
		Sudokoid(Solution solution);
        Sudokoid(char* filepath);
		
        //Solution methods
        void copySolution( const Solution solution );
        void generateRandomSolution();

		//other methods
		bool operator< (const Sudokoid &other) const;
		void Fit();
		Sudokoid Mate( Sudokoid mate, double mutationRate );
		void Print( ostream& cout);
};
