// File to test meso network forces using NVE protocol


// header files
#include "meso2D.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	bool plotCompression = 0;
	int NCELLS = 16, n1 = 24, seed = 1;
	double phi0 = 0.7, calA0 = 1.08, dispersion = 0.1, Ftol = 1e-12, Ptol = 1e-5, dt0 = 1e-2, dphi0 = 0.01;
	double boxLengthScale = 2.5;

	// name of output file
	string posf = "pos.test";
	string enf = "en.test";

	// instantiate object
	meso2D meso2Dobj(NCELLS, seed);

	// open position config file
	meso2Dobj.openPosObject(posf);
	meso2Dobj.openXtraObject(enf);

	// initialize particles are bidisperse
	meso2Dobj.initializeMesophyllCells(dispersion, calA0, phi0, Ftol, n1);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// compress to target packing pressure
	meso2Dobj.vertexJamming2D(Ftol, Ptol, dt0, dphi0, plotCompression);

	// // relax conifguration using network
	meso2Dobj.initializeMesophyllBondNetwork();
	meso2Dobj.mesoNetworkFIRE(Ftol, dt0);

	// run NVE protocol which will output configuration and energy
	double T = 1e-3;
	double ttotal = 5000.0;
	double tskip = 10.0;
	int NT = (int) floor(ttotal/dt0);
	int NPRINTSKIP = (int) floor(tskip/dt0);
	meso2Dobj.mesoNetworkNVE(T, dt0, NT, NPRINTSKIP);

	// say goodbye
	cout << "\n\n** Finished mesoNVEtest.cpp, ending. " << endl;

	return 0;
}