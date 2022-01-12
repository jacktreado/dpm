// FILE to create network of mesophyll cells in 2D
//
// * Same as mesoEnthalpyMin2D, but initial condition is read in 

// Compilation command:
// g++ -O3 --std=c++11 -I src main/meso2D/readInMesoDynMat2D.cpp src/*.cpp -o meso.o
// ./meso.o meso.input 1 1e-3 50 0.1 2 1e-4 0.5 pos.test hess.test

// header files
#include "meso2D.h"
#include <sstream>

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

// global constants
const double dispersion = 0.1; 		// polydispersity (fixed)
const double delShrink = 1e-3;		// fractional change in effective box length during extension
const double dphiPrint = 0.01;	   	// packing fractions to skip between print steps
const double boxLengthScale = 3.0; 	// neighbor list box size in units of initial l0
const double phi0 = 0.5;		   	// initial packing fraction
const double dt0 = 1e-2;		   	// initial magnitude of time step in units of MD time
const double Ftol = 1e-12; 			// force tolerance
const double dPtol = 1e-10;			// pressure change tolerance
const double phiMin = 0.3;			// minimum packing fraction in decompression algorithm
const double cL = 0.0; 				// perimeter aging (set to 0, possibly deprecated)
const double cB = 0.0; 				// bending aging (set to 0, possibly deprecated)
const double ctch = 0.5; 			// sets magnitude of \Delta U in bond-breaking
const double aL = 1.0; 				// distribution of aging to boundary (when = 1)
const double kc = 1.0; 				// interaction spring constant
const double cKb = 0; 				// change in bending energy
const int NMINSKIP = 1;				// number of frames to skip output
const int NVMAXMAG = 5; 			// scale of max number of vertices

// set parameters
const double ctcdel = 1.0;

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int seed, NVMAX;
	double kl, kb0, betaEff, da0, dl0, P0;

	// read in parameters from command line input
	string inputFile 		= argv[1];
	string kl_str 			= argv[2];
	string kb0_str 			= argv[3];
	string betaEff_str 		= argv[4];
	string da0_str 			= argv[5];
	string dl0_str 			= argv[6];
	string P0_str 			= argv[7];
	string seed_str 		= argv[8];
	string positionFile 	= argv[9];
	string hessianFile 		= argv[10];

	// using sstreams to get parameters
	stringstream klss(kl_str);
	stringstream kb0ss(kb0_str);
	stringstream betaEffss(betaEff_str);
	stringstream da0ss(da0_str);
	stringstream dl0ss(dl0_str);
	stringstream P0ss(P0_str);
	stringstream seedss(seed_str);

	// read into data
	klss >> kl;
	kb0ss >> kb0;
	betaEffss >> betaEff;
	da0ss >> da0;
	dl0ss >> dl0;
	P0ss >> P0;
	seedss >> seed;

	// instantiate object
	meso2D meso2Dobj(inputFile,seed);

	// open config & hessian file
	meso2Dobj.openPosObject(positionFile);
	meso2Dobj.openHessObject(hessianFile);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// set max # of vertices
	NVMAX = NVMAXMAG*meso2Dobj.getNVTOT();
	meso2Dobj.setNVMAX(NVMAX);

	// set aging parameters
	meso2Dobj.setbetaEff(betaEff);
	meso2Dobj.setctcdel(ctcdel);
	meso2Dobj.setctch(ctch);
	meso2Dobj.setcL(cL);
	meso2Dobj.setaL(aL);
	meso2Dobj.setcB(cB);
	meso2Dobj.setcKb(cKb);
	meso2Dobj.setkbi(kb0);
	meso2Dobj.setkl(kl);
	meso2Dobj.setkc(kc);

	// relax configuration using network + bending
	meso2Dobj.initializeMesophyllBondNetwork();
	meso2Dobj.t0ToCurrent();
	meso2Dobj.mesoFIRE(&meso2D::mesoNetworkForceUpdate, Ftol, dt0);
	// meso2Dobj.printMesoNetworkCTCS2D();
	// meso2Dobj.mesoPrintLinearResponse(&meso2D::mesoNetworkForceUpdate, Ftol, dt0);
	meso2Dobj.t0ToCurrent();

	// run stretching simulation to create network
	meso2Dobj.mesoNetworkEnthalpyMin(&meso2D::mesoNetworkForceUpdate, Ftol, dPtol, dt0, da0, dl0, P0, phiMin, NMINSKIP);

	// say goodbye
	cout << "\n** Finished readInMesoNetwork2D.cpp, ending. " << endl;

	return 0;
}