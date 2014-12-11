#include <omp.h>
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
#include "Sudokoid_8tile.h"

int GLOBAL_P;
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

	#pragma omp parallel for num_threads(GLOBAL_P) shared(matingPopulation) schedule(static, 10)
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
	#pragma omp parallel for num_threads(GLOBAL_P) shared(population) schedule(static, 10)
	for (unsigned int k = 0; k < population.size(); k++)
	{
      if(debugging)
         print_moves(population[k]);

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

returns a new population created from random fillings

Parameters:
	population_size - the size of the population to be created
******************************************************************************/
vector < Sudokoid > GenerateInitialPopulation( int population_size )
{
	vector < Sudokoid > population;
	population.resize(population_size);

	for(int i = 0; i < population_size; i++ )
	{
		population[i] = Sudokoid();
		population[i].generateRandomSolution();

      if(debugging)
         print_moves(population[i]);
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
   if(debugging)
      std::cout<< "\n---------------------Finding Best----------------\n" << std::endl;

	int lowest = INT_MAX;
	int lowestIndex = 0;

	if(debugging)
	{
		cout<<population.size();
		cout<<population[1].Fitness;
	}

	#pragma omp parallel for num_threads(GLOBAL_P) shared(population) private(lowest) schedule(static, 10)
	for(unsigned int i = 0; i < population.size(); i++)
	{
		//make sure that the fitness has been calculated
		if(population[i].Fitness == INT_MAX)
      {
			population[i].Fit();
      }

		if(population[i].Fitness < lowest)
		{
			lowest = population[i].Fitness;
			
			#pragma omp critical
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
	GLOBAL_P = 2;
	//constants, for now
	int grid_dimension = 3;
	int restart_threshold = 10;

	//set up default parameters (including genetic operators)
   //string fn = "8tiles/easy1.txt";	
    char* filename;//(char*)fn.c_str();
	int population_size = 100;
	int generations = 100;
	double selection_rate = .1;
	double mutation_rate =.05; // due to how the program is set up, must be between .0001 and 100 to work.

	//look for first few arguments
	switch(argc)
	{
		case 8: restart_threshold = atoi(argv[7]);
		case 7: mutation_rate = atof(argv[6]); //fall through to set the rest of the previous args
			if ( (mutation_rate != 0) && (mutation_rate > 100 || mutation_rate < .0001))
			{
				//because of how the mutation is calculated, any nonzero value below .0001
				//will cease to cause any effect.
				cout << "\nMutation rate must be between .0001 and 100, or 0\n";
				return 0;
			}
		case 6: selection_rate = atof(argv[5]);
		case 5: generations = atoi(argv[4]);
		case 4: population_size = atoi(argv[3]);
		case 3: filename = argv[2];
		case 2: GLOBAL_P = atoi(argv[1]);
		default:
			break;
	}
	
	//Create the first Sudokoid
	ifstream fin;

	fin.open(filename);
	if(!fin.is_open())
	{
		cout << "error opening file " << filename << endl;
		
		return 1;
	}
	
	try
	{
		readPuzzle(fin);
	}
	catch(exception& e)
	{
		cout << e.what() << endl;
	}
	//close the file
	fin.close();

	//begin output	

	cout << "Puzzle: " << filename << endl; 
	cout << "population size = " << population_size << endl;
	cout << "number of generations = " << generations << endl; 
	cout << "selection rate = " << selection_rate << endl;
	cout << "mutation rate = " << mutation_rate << endl;
	//cout << "restart threshold = " << restart_threshold << "\n\n"; 

	cout << "Initial configuration (" << grid_dimension << "x" << grid_dimension << "grid): \n\n";
	printPuzzle(cout);

	//variables for the GE
	vector < Sudokoid > population; //a vector containing the current generation of solutions
	vector < Sudokoid > champions; //a vector containing every best solution from each generation
	Sudokoid Sudoking; //the best solution so far
	Sudokoid BestSolution; //the final best solution

	champions.resize(generations);

	//generate the inital population
	int bestFit = INT_MAX;
	int generation = 0;
	population = GenerateInitialPopulation( population_size );

	int lastBest = INT_MAX; //the last best, used to find streaks.
	int tieStreak = 0; //the number of generations which have been a tie
	double startime, totalTime = 0;
	//begin the generational looping
	startime = omp_get_wtime();
	while(bestFit > 0 && generation < generations)
	{
		//create a new generation
		generation++;
		cout << "Generation " << generation << ": ";

		//select mates and breed a new generation
		vector < Sudokoid > MatingPopulation = SelectMatingPopulation(population, selection_rate);
		population = GeneratePopulation(MatingPopulation, mutation_rate, population_size);			

		//find the champion of the current generation
		Sudoking = Best(population);

		bestFit = Sudoking.Fitness;
		cout << "best score: " << bestFit << endl;

		//Add the current best (Sudoking) to the list of champions
		champions[generation - 1] = Sudoking;

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

			population = GenerateInitialPopulation( population_size );

			lastBest = INT_MAX; //the last best, used to find streaks.
			tieStreak = 0; //the number of generations which have been a tie
		}

	//select the best of all champion solutions as the best solution
	champions.resize(generation); //resize the champions so Best can tranverse it without seg faulting
	BestSolution = Best(champions);
	}
	totalTime += (omp_get_wtime() - startime);

	cout << "\nBest Solution: \n\n";
	BestSolution.Print(cout);
	
	cout << "Average Time per generation: " << totalTime/generations << endl;
		
	return 0;
}
