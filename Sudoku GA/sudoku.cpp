
/**************************************************************************//**
  Sudoku
 
 Assignemnt - Assignment 1, Genetic Algorithm
 
 Course - CSC 447 Artifical Intelligence
 
 Authors - Minda McDaniel, Derek Stotz
 
 Date - 2/24/2014
 
 Instructor - John Weiss
 
 Usage - To compile, use make
	To run, enter sudoku <command line arguments>
	At a minimum, command line arguments should be a file name, but can be
	any number in the following order:
	
		filename population_size generations selection_rate mutation_rate restart_threshold

	Defaults are as follows:
 
		population_size = 1000
		generations = 1000
		selection_rate = 0.5
		mutation_rate = .5
		restart_threshold = 100
		
Details - The purpose of this assignment was to design a sudoku solver,
	implementing a basic puzzle solver, solution evaluator, and genetic 
	algorithm to evolve a best solution.  The effectiveness of such
	a genetic algorithm depends on how well the implementation fits the
	problem and how effective the parameters used are.

	The system uses two main objects: Puzzles and Sudokoids.  A Puzzle is
	structured for solution evaluation, and the whole system begins with the
	creation of an original attempt at a solution stored in a Puzzle object.
	Sudokoids are structured for the genetic algorithm, holding 2-dimensional
	vectors of structs called SudokuCells, which are traded between mates
	during crossover.  Theses SudokuCells mutate by swaping two values
	inside of themselves.  In both Puzzles and SudokuCells, numbers in the
	puzzle are stored as characters.  Sudokoids cointain pointers to Puzzle
	counterparts, which are created out of specialized Puzzle constructors
	during the breeding proccess.

	The genetic algorithm loops as follows:

		While the generation is not the last generation,
			Increment the generation
			Select the breeding population using selection_rate
			Generate a new population by breeding the current one
			Evaluate the fitnesses of the new generation using the
				fitness function defined in the Puzzle class
			Store the most fit Sudokoid in a list of champions
		Print the best of all champions found

	The Sudokoid is set up to accept any square size of Sudoku, but
	as the Puzzle's solver is specifically for 9x9 Sudokus, this current
	implementation of the Sudokoid class only accepts 9x9 Sudoku Puzzles.
		

Recommended Usage - The more moderate difficulties can be solved quickly by:

			sudoku filename 750 500 .1 .05 10


		Interestingly enough, a mutation rate of 100% paired with a very strict
		selection rate is incredibly effective at solving the intermediate to 
		harder difficulties:
			
			sudoku filename 2000 500 .1 .1 10

		A lack of mutation compensated with a stricter selection rate and smaller
		tolerance for local minima also works suprisingly well, such as in the 
		following:

			sudoku filename 1500 500 .4 0 5
		

		All difficulties can be solved with relatively consistent speed with:

			sudoku filename 5000 500 .2 .1 10

		Larger populations and stricter selection rates give the best speed of
		success for extremely hard problems, although the generations are large
		and slow to proccess.

			sudoku filename 10000 1000 .1 .1 5



		In general:

			The population size effects the speed and effectiveness at differing rates depending
			on the difficulty of the problem.  Easier problems benefit less than harder problems
			from larger population sizes, as speedup is far more noticeable for the latter.  In
			any case, larger population sizes take much more memory.

			The maximum number of generations depends on how long the user is willing to wait
			for the algorithm to complete.  Obviously, the very hard problems should have
			smaller maximum generation counts.

			Lower selection rates seem to be more effective, and roulette wheel ranomization is
			still applied to the remaining mating population.

			The mutation rate effects populations with very high selection rates, but does not
			severely effect most of the better configurations (low selection rate, high population, low
			reset tolerance).  When it does effect the population, it helps prevent local minima.

			A reset tolerance higher than 10 is pointless, as the local maxima tend to
			restrict the algorithm quite a bit.  The resetting is very effective, however.

		
 
Issues and Bugs - running multiple sudoku proccesses on the same file simeltaneously
		will result in repeat results, as the random number generator is seeded
		on Unix time.  Runs of this program must be done at least 1 second apart.

		The local minima severely reduce the consistency of the program.
		Oftentimes extreme and fiendish difficulty problems are solved on
		the fifth generation, but sometimes it takes upwards of 40.
 
 ******************************************************************************/

#include <iostream>
#include <vector>
#include <random>
#include <time.h>
#include <fstream>
#include <exception>
#include <cstdlib>
#include <climits>
#include <algorithm>
#include <math.h>
#include "SimpleSolver.h"
#include <mpi.h>


//Deletion lists for memory managment
vector <Puzzle*> Sudokoid::DeleteList = vector <Puzzle*> ();

//Whether or not debugging messages should be displayed
bool debugging = false;

// Global variables
int ID;
int P;
int N;
/* User defined Macros */
#define BLOCK_LOW(id,p,n)	((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)	(BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n)	(BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n)	(((p)*((j)+1)-1)/(n))


/**************************************************************************************************************
GeneratePopulation

 Uses roulette wheel randomization to breed a new population of sudokoids. Returns this population in a vector.
 Assumes that the new population should be the same size as the last.
 Roulette wheel selection method inspired by
 http://stackoverflow.com/questions/10531565/how-should-roulette-wheel-selection-be-organized-for-non-sorted-population-in-g
 Random float generation method inspired by
 http://stackoverflow.com/questions/686353/c-random-float-number-generation

Parameters:
	population - the parent population's selected mating population
	mutationRate - the rate of mutation
	population_size - the size of the child population

returns - the child population
**************************************************************************************************************/
vector < Sudokoid > GeneratePopulation(vector < Sudokoid > matingPopulation, double mutationRate, int population_size)
{
	double fitnessTotal = 0.0;
	vector < double > fitnesses;
	vector < Sudokoid > children;
	fitnesses.resize(matingPopulation.size());
	children.resize(population_size);

	for(unsigned int i = 0; i < matingPopulation.size(); i++)
	{
		fitnesses[i] = 1.0 / matingPopulation[i].Fitness;
		fitnessTotal += fitnesses[i];
	}

	//set default values for the mates in case the random value lands on the last index
	Sudokoid *mateA = &matingPopulation[matingPopulation.size() - 1];
	Sudokoid *mateB = &matingPopulation[matingPopulation.size() - 1];
	bool valid; //this will flag whether or not a certain mating was valid
	int current = 0;

	//fill the child population
	while(current < population_size)
	{
		//for the first mate,
		//generate random number from 0 to fitnessTotal		
		double roulette = (fitnessTotal)*( (float)rand()/RAND_MAX);
		valid  = true;

		//subtract all values in the list until roulette passes or reaches 0
		for(unsigned int i = 0; i < matingPopulation.size(); i++)
		{
			roulette -= fitnesses[i];
			if(roulette <= 0)
			{
				//first mate has been selected!
				mateA = &matingPopulation[i];
				i = fitnesses.size(); //exit the loop gracefully
			}
		}
	
		//now obtain the same for the second mate
		roulette = (fitnessTotal)*( (float)rand()/RAND_MAX);
		for(unsigned int i = 0; i < matingPopulation.size(); i++)
		{
			roulette -= fitnesses[i];
			if(roulette <= 0)
			{
				//second mate has been selected
				mateB = &matingPopulation[i];
				//make sure there are no duplicates.
				if(mateB == mateA)
					valid = false;
				i = fitnesses.size(); //exit the loop gracefully
			}
		}

		if(valid)
		{
			// mate the two selected children			
			children[current] = mateA -> Mate(*mateB, mutationRate);	
			current++;
		}
	}
	return children;
}

/******************************************************************************
SelectMatingPopulation

removes the bottom 1 - selectionRate fraction of Sudokoids in a population, 
returning the remaining mating population. Fitness is calculated in this function,
this function does other calculations besides fitness unlike 
"SlaveSelectMatingPopulation". However, this function forks data to
other "slave" threads using MPI_Scatter. After parallel computation is complete,
MPI_Gather is used to gather the newly calculated fitness values of each
individual.


Parameters:
	population - the population of Sudokoids
	selectionRate - the selection rate

returns - the top selectionRate fraction of Sudokoids
******************************************************************************/
vector < Sudokoid > MasterSelectMatingPopulation( vector <Sudokoid> population, double selectionRate )
{
	//see if any of population needs to be fitted (should only occur on the first population)
	Puz *sol, *rsol;
	int *fitness, *rfitness;
	
	sol = (Puz*)malloc(sizeof(Puz)*N);
	rsol = (Puz*)malloc(sizeof(Puz)*N/P);
	fitness = (int*)malloc(sizeof(int)*N);
	rfitness = (int*)malloc(sizeof(int)*N/P);
	
	if(ID == 0)
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < 9; j++)
				for (int k = 0; k < 9; k++)
					for (int n = 0; n < 10; n++)
						sol[i].solution[j][k][n] = (*population[i].puzzle).puzzle.solution[j][k][n];
			fitness[i] = population[i].Fitness;
		}
	
	MPI_Scatter(sol, 9*9*10*N/P, MPI_CHAR, rsol, 9*9*10*N/P, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(fitness, N/P, MPI_INT, rfitness, N/P, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < N/P; i++)
	{
		if (rfitness[i] == INT_MAX)
			rfitness[i] = Sudokoid::Fit(rsol[i]);
	}
	MPI_Gather(rfitness, N/P, MPI_INT, fitness, N/P, MPI_INT, 0, MPI_COMM_WORLD);
	
	for(int i = 0; i < N; i++)
		population[i].Fitness = fitness[i];
	//sort a copy of the population
	vector < Sudokoid > selected = population;
	std::sort(selected.begin(), selected.end());

	//find the index at the selectionRate point
	int selectionIndex = (int)(selectionRate * selected.size());

	//change selected to only include the top Sudokoids
	selected = vector<Sudokoid> (selected.begin(), selected.begin() + selectionIndex - 1);
	return selected;

	
	
	free(sol);
	free(rsol);
	free(fitness);
	free(rfitness);
	return selected;
	
}

/******************************************************************************
SelectMatingPopulation

This function only calculates the fitness of individuals recieved from the 
data scattering from the master thread. Once done the master thread gathers
the newly calculated fitness values back to itself from this function.


Parameters:

MPI IN: sol - sub population recieved from the Master thread.
   OUT: rfitness - newly calculated fitness values of the sub population.

returns - the top selectionRate fraction of Sudokoids
******************************************************************************/
void SlaveSelectMatingPopulation( )
{
	//see if any of population needs to be fitted (should only occur on the first population)
	Puz *sol, *rsol;
	int *fitness, *rfitness;
	
	sol = (Puz*)malloc(sizeof(Puz)*N);
	rsol = (Puz*)malloc(sizeof(Puz)*N/P);
	fitness = (int*)malloc(sizeof(int)*N);
	rfitness = (int*)malloc(sizeof(int)*N/P);
	
	MPI_Scatter(sol, 9*9*10*N/P, MPI_CHAR, rsol, 9*9*10*N/P, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(fitness, N/P, MPI_INT, rfitness, N/P, MPI_INT, 0, MPI_COMM_WORLD);
	for (int i = 0; i < N/P; i++)
	{
		if (rfitness[i] == INT_MAX)
			rfitness[i] = Sudokoid::Fit(rsol[i]);
	}
	MPI_Gather(rfitness, N/P, MPI_INT, fitness, N/P, MPI_INT, 0, MPI_COMM_WORLD);

	free(sol);
	free(rsol);
	free(fitness);
	free(rfitness);	
}

/******************************************************************************
GenerateInitialPopulation

returns a new population created from random fillings of the progenitor's
empty cells.

Parameters:
	progenitor - the original best solution from the simple solver
	population_size - the size of the population to be created
******************************************************************************/
vector < Sudokoid > GenerateInitialPopulation( Sudokoid progenitor, int population_size )
{
	vector < Sudokoid > population;
	population.resize(population_size);

	for(int i = 0; i < population_size; i++ )
	{
		population[i] = Sudokoid(progenitor.Cells); //copy the progenitor
		population[i].Dimension = progenitor.Dimension;
		population[i].lockCells(); //lock all values which are not blank
		population[i].fillCells(); //fill blanks and generate new puzzles

	}

	return population;
}

/******************************************************************************
MasterBest

returns the Sudokoid with the best fitness.  Returns the first sequentially
in the case of a tie. This function transmit data to other processes for 
computaiton using MPI_Scatter and retrieves the index of the best individual
with an MPI_Reduce.

Parameter:
	population - the population to pick the best from
******************************************************************************/
Sudokoid MasterBest( vector <Sudokoid> population)
{
	int lowest = INT_MAX;
	int lowestIndex[1] = {0};
	int reducedLowest[1];

	Puz *sol, *rsol;
	int *fitness, *rfitness;

	sol = (Puz*)malloc(sizeof(Puz)*N);
	rsol = (Puz*)malloc(sizeof(Puz)*N/P);
	fitness = (int*)malloc(sizeof(int)*N);
	rfitness = (int*)malloc(sizeof(int)*N/P);
	
	if(sol == NULL || rsol == NULL || fitness == NULL || rfitness == NULL)
		MPI_Abort(MPI_COMM_WORLD, -1);
		
	int n = population.size();
	if(ID == 0)
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < 9; j++)
				for (int k = 0; k < 9; k++)
					for (int n = 0; n < 10; n++)
						sol[i].solution[j][k][n] = (*population[i].puzzle).puzzle.solution[j][k][n];
			fitness[i] = population[i].Fitness;
		}

	MPI_Scatter(sol, 9*9*10*n/P, MPI_CHAR, rsol, 9*9*10*n/P, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(fitness, n/P, MPI_INT, rfitness, n/P, MPI_INT, 0, MPI_COMM_WORLD);
	for(int i = 0; i < n/P; i++)
	{
		//make sure that the fitness has been calculated
		if(rfitness[i] == INT_MAX)
			rfitness[i] = Sudokoid::Fit(rsol[i]);

		if(rfitness[i] < lowest)
		{
			lowest = rfitness[i];
			lowestIndex[0] = i;
		}
	}
	MPI_Reduce(lowestIndex, reducedLowest, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

	free(sol);
	free(rsol);
	free(fitness);
	free(rfitness);
	return population[reducedLowest[0]];
}
/******************************************************************************
SerialBest

returns the Sudokoid with the best fitness.  Returns the first sequentially
in the case of a tie. This serial function is for the last evaluation of
the population to determine who is the best. A serial version was required
because the champions list is evaluated and the champions list length isn't
constant.

Parameter:
	population - the population to pick the best from
******************************************************************************/
Sudokoid SerialBest( vector <Sudokoid> population)
{
	int lowest = INT_MAX;
	int lowestIndex = 0;

	for(unsigned int i = 0; i < population.size(); i++)
	{
		//make sure that the fitness has been calculated
		if(population[i].Fitness == INT_MAX)
			population[i].Fit();

		if(population[i].Fitness < lowest)
		{
			lowest = population[i].Fitness;
			lowestIndex = i;
		}
	}
	return population[lowestIndex];
}

/******************************************************************************
SlaveBest

returns the Sudokoid with the best fitness.  Returns the first sequentially
in the case of a tie. This serial function is for the last evaluation of
the population to determine who is the best. This function is only executed
by slave threads/processes. It is passed no data and instead receives data
from the root thread after which the function sends the computed data
back to the root thread.

MPI IN:
	sol - the population to pick the best from
	fitness - the fitness of the population "sol"
MPI OUT:
	LowestIndex - the index of the best individual in "sol" 
******************************************************************************/
void SlaveBest()
{
	int lowest = INT_MAX;
	int lowestIndex[1] = {0};
	int reducedLowest[1];

	Puz *sol, *rsol;
	int *fitness, *rfitness;

	sol = (Puz*)malloc(sizeof(Puz)*N);
	rsol = (Puz*)malloc(sizeof(Puz)*N/P);
	fitness = (int*)malloc(sizeof(int)*N);
	rfitness = (int*)malloc(sizeof(int)*N/P);
	if(sol == NULL || rsol == NULL || fitness == NULL || rfitness == NULL)
		MPI_Abort(MPI_COMM_WORLD,-1);


	MPI_Scatter(sol, 9*9*10*N/P, MPI_CHAR, rsol, 9*9*10*N/P, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatter(fitness, N/P, MPI_INT, rfitness, N/P, MPI_INT, 0, MPI_COMM_WORLD);
	for(int i = 0; i < N/P; i++)
	{
		//make sure that the fitness has been calculated
		if(rfitness[i] == INT_MAX)
			rfitness[i] = Sudokoid::Fit(rsol[i]);

		if(rfitness[i] < lowest)
		{
			lowest = rfitness[i];
			lowestIndex[0] = i;
		}
	}
	// add offset
	lowestIndex[0] += ID*N/P;
	MPI_Reduce(lowestIndex, reducedLowest, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

	//MPI_Send

	free(sol);
	free(rsol);
	free(fitness);
	free(rfitness);
}

/******************************************************************************
main
The main function, where the top-level GE control happens and command line 
arguments are parsed and factored into the solution.

Parameters:
	argc - the number of command-line arguments
	argv - the command line arguments.

Command line arguments are as follows:
sudoku filename population_size generations selection_rate mutation_rate 

returns - an error code, but only returns 0 for now.  Errors information 
	provided via console output.
******************************************************************************/

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	MPI_Comm_rank(MPI_COMM_WORLD, &ID);
	
	//seed the random number generator
	srand ( time(NULL) );
	//ignore the first few results to avoid similarities between trials done within a short period of time.
	for (int i  = 0; i < 10; i++) { rand(); }

	//constants, for now
	int grid_dimension = 9;
	int restart_threshold = 100;
	int generation = 0;

	double time = 0.0;

	//set up default parameters (including genetic operators)
	char *filename; 
	int population_size = N = 1000;
	int generations = 1000;
	double selection_rate = 0.5;
	double mutation_rate = .05; // due to how the program is set up, must be between .0001 and 100 to work.
	//look for first few arguments
	
	int done[1] = {0};
	switch(argc)
	{
		case 7: restart_threshold = atoi(argv[6]);
		case 6: mutation_rate = atof(argv[5]); //fall through to set the rest of the previous args
			if ( (mutation_rate != 0) && (mutation_rate > 100 || mutation_rate < .0001))
			{
				//because of how the mutation is calculated, any nonzero value below .0001
				//will cease to cause any effect.
				cout << "\nMutation rate must be between .0001 and 100, or 0\n";
				return 0;
			}
		case 5: selection_rate = atof(argv[4]);
		case 4: generations = atoi(argv[3]);
		case 3: population_size = N = atoi(argv[2]);
		case 2: filename = argv[1];
			break; //done with arguments
		default: cout << "Accepts 1 to 5 arguments.\n";
			 return 0;
	}
	
	//Create the first Sudokoid
	ifstream fin;
	Sudokoid FirstPuzzleSudokoid; //The Sudokoid holding nothing but a Puzzle.  
				//Will be replaced with a Progenitor  if a GE is required for solution.
	
	if(ID == 0)
	{
		// this execution is sole for the root/master thread
		fin.open(filename);
		if(!fin.is_open())
		{
			cout << "error opening file " << filename << endl;
			return 1;
		}
		
		try
		{
			FirstPuzzleSudokoid.puzzle = new Puzzle(fin);
			//this puzzle cannot be added to DeleteList in case we need to restart
		}
		catch(exception& e)
		{
			cout << e.what() << endl;
		}
		//close the file
		fin.close();

		//begin output	

		cout << "Sudoku: " << filename << endl; 
		cout << "population size = " << population_size << endl;
		cout << "number of generations = " << generations << endl; 
		cout << "selection rate = " << selection_rate << endl;
		cout << "mutation rate = " << mutation_rate << endl;
		//cout << "restart threshold = " << restart_threshold << "\n\n"; 

		cout << "Initial configuration (" << grid_dimension << "x" << grid_dimension << "grid): \n\n";
		FirstPuzzleSudokoid.Print(cout);
		
		cout << "\n\nFilling in predetermined squares:\n\n";
		FirstPuzzleSudokoid.Fit();
		FirstPuzzleSudokoid.Print(cout);
		cout << endl;
	}
	//variables for the GE
	vector < Sudokoid > population; //a vector containing the current generation of solutions
	vector < Sudokoid > champions; //a vector containing every best solution from each generation
	Sudokoid Sudoking; //the best solution so far
	Sudokoid BestSolution; //the final best solution

	//check to see if a GE solution is even neccessary
	if(ID == 0 && FirstPuzzleSudokoid.Fitness == 0)
	{
		BestSolution = Sudokoid(FirstPuzzleSudokoid.puzzle);
	}
	else if(ID == 0)
	{
		champions.resize(generations);

		//generate the inital population
		int bestFit = INT_MAX;
		Sudokoid Progenitor = Sudokoid(FirstPuzzleSudokoid.puzzle);
		Progenitor.Dimension = sqrt ( grid_dimension ) ; // 3 x 3 cells
		population = GenerateInitialPopulation( Progenitor, N );

		int lastBest = INT_MAX; //the last best, used to find streaks.
		int tieStreak = 0; //the number of generations which have been a tie

		time = MPI_Wtime();
		//if not, begin the generational looping
		while(bestFit > 0 && generation < generations)
		{
			//create a new generation
			generation++;
			cout << "Generation " << generation << ": ";

			//select mates and breed a new generation
			vector < Sudokoid > MatingPopulation = MasterSelectMatingPopulation(population, selection_rate);
			population = GeneratePopulation(MatingPopulation, mutation_rate, N);			

			//find the champion of the current generation
			Sudoking = MasterBest(population);

			bestFit = Sudoking.Fitness;
			cout << "best score: " << bestFit << endl;

			//Add the current best (Sudoking) to the list of champions
			champions[generation - 1] = Sudoking;

			Sudokoid::deletePuzzles(); //delete puzzles to prevent memory leaks

			//test to see if there's a streak and if we need to restart.

			if(lastBest == bestFit)
			{
				//this was a tie.
				tieStreak++;
			}

			lastBest = bestFit;

			if(tieStreak >= restart_threshold)
			{
				//restart
				cout << "\nStreak of " << restart_threshold << " reached.  Restarting...\n\n";
				//generate the inital population

				Progenitor = Sudokoid(FirstPuzzleSudokoid.puzzle);
				Progenitor.Dimension = sqrt ( grid_dimension ) ; // 3 x 3 cells
				population = GenerateInitialPopulation( Progenitor, N );

				lastBest = INT_MAX; //the last best, used to find streaks.
				tieStreak = 0; //the number of generations which have been a tie
			}
			if(bestFit == 0)
				done[0] = 1;
			MPI_Bcast(done, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}

		time = MPI_Wtime() - time;

		//select the best of all champion solutions as the best solution
		champions.resize(generation); //resize the champions so Best can tranverse it without seg faulting
		BestSolution = SerialBest(champions);
		
		cout << "\nBest Solution: \n\n";
		Puzzle(&BestSolution).output(cout);
	}
	else
	{
		// slave process execute this loop
		while(generation < generations)
		{
			//create a new generation
			generation++;
			//select mates and breed a new generation
			SlaveSelectMatingPopulation();

			//find the champion of the current generation
			SlaveBest();

			Sudokoid::deletePuzzles(); //delete puzzles to prevent memory leaks
			MPI_Bcast(done, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if(done[0] != 0)
				generation = generations;
		}
		//select the best of all champion solutions as the best solution
	}


		
	//clear the heap of the beginning puzzle and any remaining puzzles
	if(ID == 0)
	{	
		if(generation > 0)
			cout << "\nAverage Time per seleting mating population and finding best individual: " 
				 << time/generation << endl;
		Sudokoid::deletePuzzles();
		delete FirstPuzzleSudokoid.puzzle;
	}
	MPI_Finalize();
	return 0;
}
