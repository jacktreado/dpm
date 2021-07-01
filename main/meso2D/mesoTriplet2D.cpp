// FILE to create mini-network with 3 mesophyll cells in 2D
// 
// * Cells are pulled to center to start
// * Contacts are formed
// * As cells are pulled away from center (no size change), shapes rigidify with aging perimeters, angles + kb
// * network is created until hmax is reached
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/meso2D/mesoTriplet2D.cpp src/*.cpp -o meso.o
// ./meso.o 24 1.01 2.0 1.0 10.0 0 0 0 1e-12 1 pos.test
// 
// 
// Parameter input list
// 1. n1: 				number of vertices on first particle
// 2. calA0: 			preferred initial shape parameter for all particles
// 3. betaEff: 			effective temperature, sets contact breaking
// 4. cL: 				perimeter aging parameter
// 5. cB: 				preferred angle aging parameter
// 6. cKb; 				bending energy aging parameter
// 7. Ftol: 			force tolerance, sets distance to each energy minimum
// 8. seed: 			seed for random number generator
// 9. positionFile: 	string of path to output file with position/configuration data
// 
// NOTE: no need to pass member function as argument, pin simulations need specific member functions


// header files
#include "meso2D.h"
#include <sstream>

// preprocessor macros
#define NDIM 2
	
// namspace
using namespace std;

// global constants
const int NCELLS 				= 3;		// always 3 cells
const double phi0 				= 0.1;		// initial packing fraction, for viz
const double dh 				= 0.0001;	// cell center step size
const double dhprint 			= 0.05;		// dh before print step
const double boxLengthScale 	= 2.5;		// neighbor list box size in units of initial l0
const double dt0 				= 1e-2;		// initial magnitude of time step in units of MD time

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int n1, seed;
	double calA0, betaEff, cL, cB, cKb, Ftol, L, hmax, kcspring;

	// read in parameters from command line input
	string n1_str 			= argv[1];
	string calA0_str 		= argv[2];
	string hmax_str 		= argv[3];
	string kcspring_str 	= argv[4];
	string betaEff_str 		= argv[5];
	string cL_str 			= argv[6];
	string cB_str 			= argv[7];
	string cKb_str 			= argv[8];
	string Ftol_str 		= argv[9];
	string seed_str 		= argv[10];
	string positionFile 	= argv[11];

	// using sstreams to get parameters
	stringstream n1ss(n1_str);
	stringstream calA0ss(calA0_str);
	stringstream hmaxss(hmax_str);
	stringstream kcspringss(kcspring_str);
	stringstream betaEffss(betaEff_str);
	stringstream cLss(cL_str);
	stringstream cBss(cB_str);
	stringstream cKbss(cKb_str);
	stringstream Ftolss(Ftol_str);
	stringstream seedss(seed_str);

	// read into data
	n1ss 			>> n1;
	calA0ss 		>> calA0;
	hmaxss			>> hmax;
	kcspringss 		>> kcspring;
	betaEffss 		>> betaEff;
	cLss 			>> cL;
	cBss 			>> cB;
	cKbss 			>> cKb;
	Ftolss 			>> Ftol;
	seedss 			>> seed;

	// instantiate object
	meso2D meso2Dobj(NCELLS, seed);

	// SET PBC -> 0, NO NEED FOR TRIPLET CELLS
	meso2Dobj.setpbc(0,false);
	meso2Dobj.setpbc(1,false);

	// open position config file
	meso2Dobj.openPosObject(positionFile);

	// initialize particles are bidisperse
	meso2Dobj.initializeMesophyllCells(0.0, calA0, phi0, Ftol, n1);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// put initial pins in box center
	L = meso2Dobj.getL(0);
	vector<double> xpin0(NDIM*NCELLS,0.0);
	fill(xpin0.begin(), xpin0.end(), 0.5*L);

	// draw pins to box center
	meso2Dobj.mesoPinFIRE(xpin0, Ftol, dt0, kcspring);

	// set aging parameters
	meso2Dobj.setbetaEff(betaEff);
	meso2Dobj.setcL(cL);
	meso2Dobj.setcB(cB);
	meso2Dobj.setcKb(cKb);

	// initialize adhesive network contacts
	meso2Dobj.initializeMesophyllBondNetwork();

	// run stretching simulation to create network
	meso2Dobj.mesoPinExtension(Ftol, dt0, hmax, dh, dhprint, 2.0*kcspring);


	// say goodbye
	cout << "\n** Finished mesoTriplet2D.cpp, ending. " << endl;

	return 0;
}