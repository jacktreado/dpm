// File to jam configuration of DPM particles to target pressure
// 
// Will create bidisperse DPM particles, set constants,
// place particle centers, relax shapes + positions, compress to target pressure Ptol (jamming is limit of Ptol -> 0), print configuration
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/jamtest.cpp src/*.cpp -o test.o
// 
// Run command:
// ./test.o


// header files
#include "dpm.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int NCELLS = 12, nsmall = 24, seed = 1;
	double phi0 = 0.2, calA0 = 1.2, smallfrac = 0.5, sizefrac = 1.4, Ftol = 1e-12, Ptol = 1e-8, dt0 = 1e-2;
	double ka = 1.0, kl = 1.0, kb = 0.01, kc = 1.0, thA = 12.0, thK = 3.0, boxLengthScale = 2.5, l1 = 0.0, l2 = 0.0;

	// options for attraction
	bool useAttraction = 0;

	// pointer to dpm member function (should pt to null)
	dpmMemFn forceUpdate = nullptr;

	// name of output file
	string posf = "pos.test";

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

	// compress to target packing fraction
	double phi0Target = 1.0, dphi0 = 0.005;
	bool plotCompression = true;
	configobj2D.vertexJamming2D(forceUpdate,Ftol,Ptol,dt0,dphi0,plotCompression);

	// say goodbye
	cout << "\n\n** Finished jamtest.cpp, ending. " << endl;

	return 0;
}