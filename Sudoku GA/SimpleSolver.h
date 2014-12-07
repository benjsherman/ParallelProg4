#include<iostream>
#include"puzzle.h"


using namespace std;

class SimpleSolver
{
	public:
		static int ssolve(char sol[9][9][10]);
		static int getSolutionScore(char sol[9][9][10]);

	private:
		static void markSqrPossiblities(int x, int y, char sol[9][9][10]);
		static void hiddenSingleSquare(int x, int y, Puzzle &p);
		static void resetPossiblities(char sol[9][9][10]);
		static void mrkRowColmPossiblities(char sol[9][9][10]);
		static void findNakedSingles(char sol[9][9][10]);
		static void findHiddenSingleRowColn(char sol[9][9][10]);
		static int squareDuplicates(int x, int y, char sol[9][9][10]);
		static int rowDuplicates(char sol[9][9][10]);
		static int countBlanks(char sol[9][9][10]);
		static int colDuplicates(char sol[9][9][10]);
		static void hiddenSingleSquare(int x, int y, char sol[9][9][10]);
};
