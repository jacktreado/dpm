// File to model crawling of multiple cells on plate
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/crawling_test.cpp src/*.cpp -o test.o
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
	int NCELLS = 16, tNV = 20, seed = 1;
	int NT = 1e6, NPRINTSKIP = 5e3;
	double tCalA0 = 1.15, tDisp = 0.15;
	double phi0 = 0.1, dt0 = 2e-2, Ftol = 1e-12, kwell = 1e-2;
	double ka = 1.0, kl = 1.0, kb = 0.0, kc = 1.0, gamtt = -0.5, boxLengthScale = 1.5, l1 = 0, l2 = 0;
	double v0 = 0.1, Dr0 = 0.5, Ds = 0.3;
	tumor2DMemFn crawlingForceUpdate = nullptr;

	// name of output file
	string posf = "pos.test";

	// instantiate object
	tumor2D tumor2Dobj(NCELLS, seed);

	// set periodic boundaries
	tumor2Dobj.setpbc(0,1);
	tumor2Dobj.setpbc(1,1);

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
	tumor2Dobj.setl1(l1);
	tumor2Dobj.setl2(l2);

	// initialize tumor cells
	tumor2Dobj.gaussian2D(tDisp, tCalA0, tNV);

	// initialize particle positions
	tumor2Dobj.initializeTumorMonolayerPositions(phi0, Ftol, kwell);

	// initialize neighbor linked list
	tumor2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// force updates for crawling
	crawlingForceUpdate = &tumor2D::stickyTumorForceUpdate;

	// relax positions using tumor FIRE
	tumor2Dobj.tumorFIRE(crawlingForceUpdate, Ftol, dt0);

	// -- crawling simming
	tumor2Dobj.setgamtt(gamtt);
	tumor2Dobj.setdt(dt0);
	tumor2Dobj.crawling(crawlingForceUpdate,NT,NPRINTSKIP);

	// say goodbye
	cout << "\n\n** Finished invasion.cpp, ending. " << endl;

	return 0;
}