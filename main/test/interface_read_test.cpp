// File to test read in to interface
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/interface_read_test.cpp src/*.cpp -o test.o
// 
// Run command:
// ./test.o

// header files
#include "tumor2D.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int seed = 1;
	double ka = 1.0, kl = 1.0, kb = 0.0, kc = 1.0, boxLengthScale = 2.5, l1 = 0.01, l2 = 0.05;
	double v0 = 0.03, Dr0 = 0.1, Ds = 0.2, ecmbreak = 0.5, kecm = 0.1;

	// name of output file
	string posf = "pos_check.test";
	string inputf = "input.test";

	// instantiate object
	tumor2D tumor2Dobj(inputf,seed);

	// open position config file
	tumor2Dobj.openPosObject(posf);

	// set spring constants
	tumor2Dobj.setka(ka);
	tumor2Dobj.setkl(kl);
	tumor2Dobj.setkb(kb);
	tumor2Dobj.setkc(kc);

	// activity parameters
	tumor2Dobj.setv0(v0);
	tumor2Dobj.setDr0(Dr0);
	tumor2Dobj.setDs(Ds);
	tumor2Dobj.setkecm(kecm);
	tumor2Dobj.setecmbreak(ecmbreak);

	// initialize neighbor linked list
	tumor2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// print configuration back to pos.test
	tumor2Dobj.printTumorInterface(0.0);

	// say goodbye
	cout << "\n\n** Finished interface_read_test.cpp, ending. " << endl;

	return 0;
}