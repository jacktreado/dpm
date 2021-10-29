// FILE to create network of mesophyll cells in 2D
//
// * Same as mesoEnthalpyMin2D, but initial condition is read in 

// Compilation command:
// g++ -O3 --std=c++11 -I src main/meso2D/readInMesoNetwork2D.cpp src/*.cpp -o meso.o
// ./meso.o meso.input 1e-2 100 1e-3 1.5 1e-6 0.5 1 1 1 pos.test
//
//
// Parameter input list
// 1. inputFile			file to 
// 2. kb0: 				initial amount of bending energy
// 3. betaEff: 			inverse bond breaking temperature
// 4. ctcdel: 			cda (1) or not (0)
// 5. cL: 				perimeter aging parameter
// 6. aL: 				distribution of aging to either contacts (0) or void (1)
// 7. cB: 				preferred angle aging parameter
// 8. seed: 			seed for random number generator
// 9. positionFile: 	string of path to output file with position/configuration data

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
const double boxLengthScale = 2.5; 	// neighbor list box size in units of initial l0
const double phi0 = 0.5;		   	// initial packing fraction
const double dt0 = 5e-3;		   	// initial magnitude of time step in units of MD time
const double Ftol = 1e-10; 			// force tolerance
const double dPtol = 1e-10;			// pressure change tolerance
const double phiMin = 0.3;			// minimum packing fraction in decompression algorithm
const double kl = 0.5; 				// perimeter spring constant
const double aL = 1.0; 				// distribution of aging to boundary (when = 1)
const double kc = 0.5; 				// interaction spring constant
const double cKb = 0; 				// change in bending energy
const int NMINSKIP = 10;			// number of frames to skip output
const int NVMAXMAG = 5; 			// scale of max number of vertices

// set parameters
const double ctcdel = 1.0;

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int seed, NVMAX;
	double kb0, betaEff, da0, dl0, P0, ctch, cL, cB;

	// read in parameters from command line input
	string inputFile 		= argv[1];
	string kb0_str 			= argv[2];
	string betaEff_str 		= argv[3];
	string da0_str 			= argv[4];
	string dl0_str 			= argv[5];
	string P0_str 			= argv[6];
	string ctch_str 		= argv[7];
	string cL_str 			= argv[8];
	string cB_str 			= argv[9];
	string seed_str 		= argv[10];
	string positionFile 	= argv[11];

	// using sstreams to get parameters
	stringstream kb0ss(kb0_str);
	stringstream betaEffss(betaEff_str);
	stringstream da0ss(da0_str);
	stringstream dl0ss(dl0_str);
	stringstream P0ss(P0_str);
	stringstream ctchss(ctch_str);
	stringstream cLss(cL_str);
	stringstream cBss(cB_str);
	stringstream seedss(seed_str);

	// read into data
	kb0ss >> kb0;
	betaEffss >> betaEff;
	da0ss >> da0;
	dl0ss >> dl0;
	P0ss >> P0;
	ctchss >> ctch;
	cLss >> cL;
	cBss >> cB;
	seedss >> seed;

	// instantiate object
	meso2D meso2Dobj(inputFile,seed);

	// open position config file
	meso2Dobj.openPosObject(positionFile);

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
	meso2Dobj.t0ToCurrent();
	meso2Dobj.printMesoNetworkCTCS2D();

	// run stretching simulation to create network
	meso2Dobj.mesoNetworkEnthalpyMin(&meso2D::mesoNetworkForceUpdate, Ftol, dPtol, dt0, da0, dl0, P0, phiMin, NMINSKIP);

	// say goodbye
	cout << "\n** Finished readInMesoNetwork2D.cpp, ending. " << endl;

	return 0;
}