// File to test initialization + FIRE relaxation of DPM particles
// 
// 	Will create bidisperse DPM particles, set constants,
// 	place particle centers, relax shapes + positions, print configuration
//  g++ -O3 -std=c++11 -I src main/test/firetest.cpp src/*.cpp -o firetest.o 



// header files
#include "dpm.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int NCELLS = 16, nsmall = 24, seed = 1;
	double phi0 = 0.7, calA0 = 1.2, smallfrac = 0.5, sizefrac = 1.4, Ftol = 1e-12, dt0 = 1e-2;
	double dispersion = 0.1;
	double ka = 1.0, kl = 1.0, kb = 0.01, kc = 1.0, thA = 10.0, thK = 3.0, boxLengthScale = 2.5;

	// name of output file
	string posf = "pos.test";

	// instantiate object
	dpm configobj2D(NCELLS, seed);

	// open position config file
	configobj2D.openPosObject(posf);

	// set spring constants
	configobj2D.setka(ka);
	configobj2D.setkl(kl);
	configobj2D.setkb(kb);
	configobj2D.setkc(kc);

	// initialize particles are bidisperse
	// configobj2D.bidisperse2D(calA0, nsmall, smallfrac, sizefrac);
	configobj2D.gaussian2D(dispersion, calA0, nsmall);

	// set preferred angle to sinusoidal
	configobj2D.sinusoidalPreferredAngle(thA, thK);

	// initialize particle positions
	configobj2D.initializePositions2D(phi0, Ftol);

	// initialize neighbor linked list
	configobj2D.initializeNeighborLinkedList2D(boxLengthScale);

	// use FIRE to relax configuration interacting via bumps to energy minimum
	configobj2D.vertexFIRE2D(&dpm::forceUpdate, Ftol, dt0);

	// print config
	configobj2D.printConfiguration2D();

	// say goodbye
	cout << "\n\n** Finished firetest.cpp, ending. " << endl;

	return 0;
}