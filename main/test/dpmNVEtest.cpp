// File to test dpm forces using NVE protocol
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/dpmNVEtest.cpp src/*.cpp -o test.o

// header files
#include "dpm.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int NCELLS = 12, nsmall = 32, seed = 1;
	double phi0 = 0.6, calA0 = 1.2, smallfrac = 0.5, sizefrac = 1.4, Ftol = 1e-12, Ptol = 1e-8, dt0 = 1e-2;
	double ka = 1.0, kl = 1.0, kb = 0.1, kc = 1.0, thA = 12.0, thK = 3.0, boxLengthScale = 2.5, l1 = 0.0, l2 = 0.0;

	// options for attraction
	bool useAttraction = 1;

	// pointer to dpm member function (should pt to null)
	dpmMemFn forceUpdate = nullptr;

	// name of output file
	string posf = "pos.test";
	string enf = "en.test";

	// open energy file in main
	ofstream enout(enf.c_str());
	if (!enout.is_open()){
		cerr << "\t** ERROR: Energy file " << enf << " could not open, ending here." << endl;
		return 1;
	}

	// instantiate object
	dpm configobj2D(NCELLS, NDIM, seed);

	// open position config file
	configobj2D.openPosObject(posf);

	// set spring constants
	configobj2D.setka(ka);
	configobj2D.setkl(kl);
	configobj2D.setkb(kb);
	configobj2D.setkc(kc);

	if (useAttraction) {
		// attractive force update
		forceUpdate = &dpm::attractiveForceUpdate;
		l1 = 0.01;
		l2 = 0.02;
		configobj2D.setl1(l1);
		configobj2D.setl2(l2);
	}
	else
		forceUpdate = &dpm::repulsiveForceUpdate;

	// initialize particles are bidisperse
	// configobj2D.bidisperse2D(calA0, nsmall, smallfrac, sizefrac);
	configobj2D.gaussian2D(0.0, calA0, nsmall);

	// set preferred angle to sinusoidal
	configobj2D.sinusoidalPreferredAngle(thA, thK);

	// initialize particle positions
	configobj2D.initializePositions2D(phi0, Ftol);

	// initialize neighbor linked list
	configobj2D.initializeNeighborLinkedList2D(boxLengthScale);

	// run NVE protocol which will output configuration and energy
	double T = 1e-4;
	double ttotal = 1000.0;
	double tskip = 10.0;
	int NT = (int) floor(ttotal/dt0);
	int NPRINTSKIP = (int) floor(tskip/dt0);
	configobj2D.vertexNVE2D(enout, forceUpdate, T, dt0, NT, NPRINTSKIP);

	// say goodbye
	cout << "\n\n** Finished dpmNVEtest.cpp, ending. " << endl;

	// close file
	enout.close();

	return 0;
}