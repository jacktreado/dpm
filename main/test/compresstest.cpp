// File to test initialization + FIRE relaxation + compression of DPM particles
// 
// 	Will create bidisperse DPM particles, set constants,
// 	place particle centers, relax shapes + positions, compress to target phi, print configuration
//  g++ -O3 --std=c++11 -I src main/test/compresstest.cpp src/*.cpp -o compresstest.o 



// header files
#include "dpm.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int NCELLS = 24, nsmall = 24, seed = 1;
	double phi0 = 0.3, calA0 = 1.17, smallfrac = 0.5, sizefrac = 1.4, disp = 0.1, Ftol = 1e-12, dPtol = 1e-12, dt0 = 1e-2;
	double ka = 1.0, kl = 1.0, kb = 0.0, kc = 1.0, thA = 0.0, thK = 0.0, boxLengthScale = 3.5;
	double P0 = 1e-6, dphi0 = 0.01;
	bool plotCompression = 1;

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
	configobj2D.bidisperse2D(calA0, nsmall, smallfrac, sizefrac);
	// configobj2D.gaussian2D(disp, calA0, nsmall);

	// set preferred angle to sinusoidal
	configobj2D.sinusoidalPreferredAngle(thA, thK);

	// initialize particle positions
	configobj2D.initializePositions2D(phi0, Ftol);

	// initialize neighbor linked list
	configobj2D.initializeNeighborLinkedList2D(boxLengthScale);

	// compress to target packing fraction
	// configobj2D.vertexCompress2Target2D(&dpm::repulsiveForceUpdate,Ftol,dt0,phi0Target,dphi0);
	configobj2D.vertexEnthalpyMin(&dpm::repulsiveForceUpdate,Ftol,dPtol,P0,dt0,plotCompression);

	// say goodbye
	cout << "\n\n** Finished compresstest.cpp, ending. " << endl;

	return 0;
}