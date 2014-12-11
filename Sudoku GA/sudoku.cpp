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
returning the remaining mating population.

Parameters:
	population - the population of Sudokoids
	selectionRate - the selection rate

returns - the top selectionRate fraction of Sudokoids
******************************************************************************/
vector < Sudokoid > SelectMatingPopulation( vector <Sudokoid> population, double selectionRate )
{
	//see if any of population needs to be fitted (should only occur on the first population)
	for (unsigned int k = 0; k < population.size(); k++)
	{
		if(population[k].Fitness == INT_MAX)
			population[k].Fit();
	}

	//sort a copy of the population
	vector < Sudokoid > selected = population;
	std::sort(selected.begin(), selected.end());

	//find the index at the selectionRate point
	int selectionIndex = (int)(selectionRate * selected.size());

	//change selected to only include the top Sudokoids
	selected = vector<Sudokoid> (selected.begin(), selected.begin() + selectionIndex - 1);
	return selected;
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
Best

returns the Sudokoid with the best fitness.  Returns the first sequentially
in the case of a tie.

Parameter:
	population - the population to pick the best from
******************************************************************************/
Sudokoid Best( vector <Sudokoid> population)
{
	int lowest = INT_MAX;
	int lowestIndex = 0;

	if(debugging)
	{
		cout<<population.size();
		cout<<population[1].Fitness;
		return Sudokoid();
	}

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
	//seed the random number generator
	srand ( time(NULL) );
	//ignore the first few results to avoid similarities between trials done within a short period of time.
	for (int i  = 0; i < 10; i++) { rand(); }

	//constants, for now
	int grid_dimension = 9;
	int restart_threshold = 100;

	//set up default parameters (including genetic operators)
	char *filename; 
	int population_size = 1000;
	int generations = 1000;
	double selection_rate = 0.5;
	double mutation_rate = .05; // due to how the program is set up, must be between .0001 and 100 to work.

	// MPI variables 
	int my_rank, processes;
	MPI_Comm comm = MPI_COMM_WORLD;
	int i = 0;	
	
	//look for first few arguments
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
		case 3: population_size = atoi(argv[2]);
		case 2: filename = argv[1];
			break; //done with arguments
		default: cout << "Accepts 1 to 5 arguments.\n";
			 return 0;
	}
	
	//Create the first Sudokoid
	ifstream fin;
	Sudokoid FirstPuzzleSudokoid; //The Sudokoid holding nothing but a Puzzle.  
				//Will be replaced with a Progenitor  if a GE is required for solution.
	
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

	//variables for the GE
	vector < Sudokoid > population; //a vector containing the current generation of solutions
	vector < Sudokoid > champions; //a vector containing every best solution from each generation
	Sudokoid Sudoking; //the best solution so far
	Sudokoid BestSolution; //the final best solution

	//check to see if a GE solution is even neccessary
	if(FirstPuzzleSudokoid.Fitness == 0)
	{
		BestSolution = Sudokoid(FirstPuzzleSudokoid.puzzle);
	}
	else
	{
		
		champions.resize(generations);

		//generate the inital population
		int bestFit = INT_MAX;
		int generation = 0;
		Sudokoid Progenitor = Sudokoid(FirstPuzzleSudokoid.puzzle);
		Progenitor.Dimension = sqrt ( grid_dimension ) ; // 3 x 3 cells
		population = GenerateInitialPopulation( Progenitor, population_size );

		int lastBest = INT_MAX; //the last best, used to find streaks.
		int tieStreak = 0; //the number of generations which have been a tie
		int cont = 1;
		MPI_Status status;


		int recv = 0;
							
		MPI_Init(NULL, NULL);
		
		MPI_Comm_size (comm, &processes);
		MPI_Comm_rank (comm, &my_rank);	
		
		if (my_rank == 0)
		{	
			// Champs received from each thread are stored here
			while (!recv)
				MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &recv, &status);
			
			// gets to here when a 0 is sent, may not be receving or sending
			int sender = status.MPI_SOURCE;

			MPI_Recv (&cont, 1, MPI_INT, sender, MPI_ANY_TAG, comm, &status);
			for (i = 1; i < processes; i++)  // may not be sending or others may not be receiveing
				MPI_Send (&recv, 1, MPI_INT, i, 0, comm); 
		}	
		//if not, begin the generational looping
		else if (my_rank != 0)
		{
			while(bestFit > 0 && generation < generations)
			{
				//create a new generation
				generation++;

				cout << "MY_RANK: " << my_rank << " Generation " << generation << ": ";

				//select mates and breed a new generation
				vector < Sudokoid > MatingPopulation = SelectMatingPopulation(population, selection_rate);
				population = GeneratePopulation(MatingPopulation, mutation_rate, population_size);			

				//find the champion of the current generation
				Sudoking = Best(population);

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
					population = GenerateInitialPopulation( Progenitor, population_size );

					lastBest = INT_MAX; //the last best, used to find streaks.
					tieStreak = 0; //the number of generations which have been a tie
				}
				MPI_Iprobe (0, MPI_ANY_TAG, comm, &recv, &status);
				if (recv)
					break;	
				
			}

			if (bestFit == 0 )
			{
				cont = 0;			
				MPI_Send(&cont, 1, MPI_INT, 0, 10, comm);
			}
			// else gather all champions and find best, store in best solution 	 
		}
			
		//select the best of all champion solutions as the best solution
		if (my_rank != 0)
		{
			champions.resize(generation); //resize the champions so Best can tranverse it without seg faulting
			BestSolution = Best(champions);
		}
	}

	// MPI Reduce the fit values of the best solutions, get the thread with lowest fit, print the solution from the rank with the lowest reduce value
	int localResult[2];

	if (my_rank == 0)
		localResult[0] = INT_MAX;
	else 
		localResult[0] = BestSolution.Fitness;
	localResult[1] = my_rank;

	int globalResult[2] = { 0, 0 };

	MPI_Allreduce(localResult, globalResult, 1, MPI_2INT, MPI_MINLOC, comm);
	
	if (my_rank == globalResult[1])
	{
		cout << "\nBest Solution found by process " << my_rank << ": \n\n";
		Puzzle(&BestSolution).output(cout);
	}

	MPI_Finalize();

	//clear the heap of the beginning puzzle and any remaining puzzles
	Sudokoid::deletePuzzles();
	delete FirstPuzzleSudokoid.puzzle;

	return 0;
}
