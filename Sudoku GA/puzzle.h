#include <iostream>
#include "Sudokoid.h"

using namespace std;


/*
Model Class for Puzzle
*/
class Puzzle
{
	public:
		Puzzle(istream & _in);
		Puzzle(Sudokoid* const sudokoid);
		Puzzle(Puzzle* puzzle);
		void output(ostream & _out);
		char solution[9][9][10]; //third dimension used to mark possible answers
		int getSolutionScore();
	private:
		int solScore; //number of blank or duplicate values (0 is perfect score)
		int countBlanks();
		int rowDuplicates();
		int colDuplicates();
		int squareDuplicates(int x, int y);
};
