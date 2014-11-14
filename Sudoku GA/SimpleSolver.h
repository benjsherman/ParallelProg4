#include<iostream>
#include"puzzle.h"

using namespace std;

class SimpleSolver
{
	public:
		static void ssolve(Puzzle &p);
	private:
		static void markSqrPossiblities(int x, int y, Puzzle &p);
		static void hiddenSingleSquare(int x, int y, Puzzle &p);
		static void resetPossiblities(Puzzle &p);
		static void mrkRowColmPossiblities(Puzzle &p);
		static void findNakedSingles(Puzzle &p);
		static void findHiddenSingleRowColn(Puzzle &p);
};
