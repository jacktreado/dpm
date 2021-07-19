// File to test tumor invasion
// 
// Will create tumor/adipocyte DPM particles, set constants,
// place particle centers, relax shapes + positions, compress to target pressure Ptol (jamming is limit of Ptol -> 0),
// invade for a bit, print data to file
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/invasiontest.cpp src/*.cpp -o test.o
// 
// Run command:
// ./test.o
// 
// TO DO (07/15)
// -- Check areaRatio, seems to be too large, slowing down compression algorithm?

// header files
#include "tumor2D.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int NCELLS, aN = 6, tN, aNV = 24, tNV = 24, seed = 1;
	int NT = 1e5, NPRINTSKIP = 1e3;
	double areaRatio = 4.0, aCalA0 = 1.01, tCalA0 = 1.20, aDisp = 0.1, tDisp = 0.1, prt = 1.0/3.0;
	double phi0 = 0.3, dphi0 = 0.01, Ftol = 1e-12, Ptol = 1e-5, dt0 = 2e-2;
	double ka = 1.0, kl = 1.0, kb = 0.0, kc = 1.0, boxLengthScale = 2.5, l1 = 0.01, l2 = 0.05;
	double v0 = 0.03, Dr0 = 0.1, Ds = 0.2, dDr = 0.1, dPsi = 0.1, Drmin = 1e-4, ecmbreak = 0.5, kecm = 0.1;
	bool plotcompression = 0;
	tumor2DMemFn invasionForceUpdate = nullptr;

	// determine number of tumor cells based on areaRatio and prt
	tN = round(aN * areaRatio * (prt/(1 - prt)));
	NCELLS = tN + aN;

	// name of output file
	string posf = "pos.test";

	// instantiate object
	tumor2D tumor2Dobj(NCELLS, tN, seed);

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

	// initialize adipocyte and tumor cells
	tumor2Dobj.initializeTumorInterface(aCalA0, tCalA0, aDisp, tDisp, areaRatio, aNV, tNV);

	// initialize particle positions
	tumor2Dobj.initializeTumorInterfacePositions(phi0, Ftol, prt);

	// initialize neighbor linked list
	tumor2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// compress to jamming at particular pressure

	// force updates for jamming + invasion
	invasionForceUpdate = &tumor2D::stickyTumorInterfaceForceUpdate;

	// -- COMPRESSION attraction parameters
	tumor2Dobj.tumorCompression(Ftol,Ptol,dt0,dphi0);
	return 0;

	// invasion simulation

	// -- INVASION attraction parameters
	tumor2Dobj.setl1(l1);
	tumor2Dobj.setl2(l2);
	tumor2Dobj.setdt(dt0);
	tumor2Dobj.invasion(invasionForceUpdate,dDr,dPsi,Drmin,NT,NPRINTSKIP);

	// say goodbye
	cout << "\n\n** Finished invasion.cpp, ending. " << endl;

	return 0;
}