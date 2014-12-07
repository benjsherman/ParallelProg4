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
		Puzzle(Puzzle* puz);
		void output(ostream & _out);
		Puz puzzle; //third dimension used to mark possible answers
	private:
};
