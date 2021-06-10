// File to test initialization of DPM particles
// 
// 	Will create bidisperse DPM particles, set constants,
// 	place particle centers, print vertex locations




// header files
#include "dpm.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int NCELLS = 16, nsmall = 24, seed = 1;
	double phi0 = 0.7, calA0 = 1.02, smallfrac = 0.5, sizefrac = 1.4, Ftol = 1e-12;
	double ka = 1.0, kl = 1.0, kb = 1.0, kc = 1.0, thA = 1.0, thK = 3.0, boxLengthScale = 2.5;

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
	configobj2D.bidisperse2D(calA0, nsmall, smallfrac, sizefrac);

	// set preferred angle to sinusoidal
	configobj2D.sinusoidalPreferredAngle(thA, thK);

	// initialize particle positions
	configobj2D.initializePositions2D(phi0, Ftol);

	// initialize neighbor linked list
	configobj2D.initializeNeighborLinkedList2D(boxLengthScale);

	// compute shape forces
	configobj2D.updateForces();

	// print config
	configobj2D.printConfiguration2D();

	// print out all net forces, to check
	int NVTOT = configobj2D.getNVTOT();
	int gi;
	double fxsum, fysum;
	for (int ci=0; ci<NCELLS; ci++){
		fxsum = 0.0, fysum = 0.0;
		for (int vi=0; vi<configobj2D.getNV(ci); vi++){
			gi = configobj2D.gindex(ci,vi);
			cout << "gi = " << gi << "  (" << NVTOT << ");  ";
			cout << "ci = " << ci << "  (" << NCELLS << "),  " << "vi = " << vi << "  (" << configobj2D.getNV(ci) << "),  ";
			cout << "F: " << configobj2D.getF(gi,0) << ",  " << configobj2D.getF(gi,1) << endl;
			fxsum += configobj2D.getF(gi,0);
			fysum += configobj2D.getF(gi,1);
		}
		cout << "\t** fxsum = " << fxsum << ";  fysum = " << fysum << endl;
	}

	// print contact matrix to check
	configobj2D.printContactMatrix();

	// say goodbye
	cout << "\n\n** Finished init.cpp, ending. " << endl;

	return 0;
}