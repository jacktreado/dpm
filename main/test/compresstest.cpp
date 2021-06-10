// File to test initialization + FIRE relaxation + compression of DPM particles
// 
// 	Will create bidisperse DPM particles, set constants,
// 	place particle centers, relax shapes + positions, compress to target phi, print configuration




// header files
#include "dpm.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int NCELLS = 16, nsmall = 28, seed = 1;
	double phi0 = 0.3, calA0 = 1.4, smallfrac = 0.5, sizefrac = 1.4, Ftol = 1e-12, dt0 = 1e-2;
	double ka = 1.0, kl = 1.0, kb = 0.01, kc = 1.0, thA = 12.0, thK = 3.0, boxLengthScale = 2.5;

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
	configobj2D.compress2Target(Ftol,dt0,phi0Target,dphi0);

	// say goodbye
	cout << "\n\n** Finished compresstest.cpp, ending. " << endl;

	return 0;
}