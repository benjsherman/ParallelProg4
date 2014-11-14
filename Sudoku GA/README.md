Program: Sudoku Solver

Last Modified: Feb 24, 2014

Compiling: No special compilation instructions; simply run make.

Usage:
	sudoku filename population_size generations selection_rate mutation_rate
	
	The only required argument is filename. The file must contain a 9 x 9 
	sudoku puzzle. We did not attempt the extra credit for other sizes.
	
	The defaults for the optional arguments are as follows:
		population_size = 1000
		generations = 1000
		selection_rate = 0.5
		mutation_rate = .05
	
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
	during crossover.  Theses SudokuCells mutate by swapping two values
	inside of themselves.  In both Puzzles and SudokuCells, numbers in the
	puzzle are stored as characters.  Sudokoids contain pointers to Puzzle
	counterparts, which are created out of specialized Puzzle constructors
	during the breeding process.

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

			Lower selection rates seem to be more effective, and roulette wheel randomization is
			still applied to the remaining mating population.

			The mutation rate effects populations with very high selection rates, but does not
			severely effect most of the better configurations (low selection rate, high population, low
			reset tolerance).  When it does effect the population, it helps prevent local minima.

			A reset tolerance higher than 10 is pointless, as the local maxima tend to
			restrict the algorithm quite a bit.  The resetting is very effective, however.

		
 
Issues and Bugs - running multiple sudoku processes on the same file simultaneously
		will result in repeat results, as the random number generator is seeded
		on Unix time.  Runs of this program must be done at least 1 second apart.

		The local minima severely reduce the consistency of the program.
		Often times extreme and fiendish difficulty problems are solved on
		the fifth generation, but sometimes it takes upwards of 40.