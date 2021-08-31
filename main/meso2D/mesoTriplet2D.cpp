// FILE to create mini-network with 3 mesophyll cells in 2D
// 
// * Cells are pulled to center to start
// * Contacts are formed
// * As cells are pulled away from center (no size change), shapes rigidify with aging perimeters, angles + kb
// * network is created until hmax is reached
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/meso2D/mesoTriplet2D.cpp src/*.cpp -o meso.o
// ./meso.o 24 1.01 1e-4 0.1 10.0 0.005 1.0 0 0 1 pos.test
// 
// 
// Parameter input list
// 1. n1: 				number of vertices on first particle
// 2. calA0: 			preferred initial shape parameter for all particles
// 3. dh 				step size
// 4. kcspring 			pin spring constant
// 5. betaEff: 			effective temperature, sets contact breaking
// 6. cL: 				perimeter aging parameter
// 7. aL: 				distribution of aging parameter to contact (0) vs void (1)
// 8. cB: 				preferred angle aging parameter
// 9. cKb; 				bending energy aging parameter
// 10. seed: 			seed for random number generator
// 11. positionFile: 	string of path to output file with position/configuration data
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
const double hmax 				= 2.5;		// max step length
const double dhprint 			= 0.05;		// dh before print step
const double boxLengthScale 	= 2.5;		// neighbor list box size in units of initial l0
const double dt0 				= 1e-2;		// initial magnitude of time step in units of MD time
const double Ftol 				= 1e-12; 	// force tolerance
const double kb0 				= 1e-4; 	// initial bending energy

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int n1, seed;
	double calA0, betaEff, cL, aL, cB, cKb, L, dh, kcspring;

	// read in parameters from command line input
	string n1_str 			= argv[1];
	string calA0_str 		= argv[2];
	string dh_str 			= argv[3];
	string kcspring_str 	= argv[4];
	string betaEff_str 		= argv[5];
	string cL_str 			= argv[6];
	string aL_str 			= argv[7];
	string cB_str 			= argv[8];
	string cKb_str 			= argv[9];
	string seed_str 		= argv[10];
	string positionFile 	= argv[11];

	// using sstreams to get parameters
	stringstream n1ss(n1_str);
	stringstream calA0ss(calA0_str);
	stringstream dhss(dh_str);
	stringstream kcspringss(kcspring_str);
	stringstream betaEffss(betaEff_str);
	stringstream cLss(cL_str);
	stringstream aLss(aL_str);
	stringstream cBss(cB_str);
	stringstream cKbss(cKb_str);
	stringstream seedss(seed_str);

	// read into data
	n1ss 			>> n1;
	calA0ss 		>> calA0;
	dhss			>> dh;
	kcspringss 		>> kcspring;
	betaEffss 		>> betaEff;
	cLss 			>> cL;
	aLss 			>> aL;
	cBss 			>> cB;
	cKbss 			>> cKb;
	seedss 			>> seed;

	// instantiate object
	meso2D meso2Dobj(NCELLS, seed);

	// SET PBC -> 0, NO NEED FOR TRIPLET CELLS
	meso2Dobj.setpbc(0,false);
	meso2Dobj.setpbc(1,false);

	// open position config file
	meso2Dobj.openPosObject(positionFile);

	// initialize particles
	meso2Dobj.initializeMesophyllCells(0.0, calA0, phi0, Ftol, n1);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// put initial pins in box center
	L = meso2Dobj.getL(0);
	vector<double> xpin0(NDIM*NCELLS,0.0);
	double th = 0.0;
	for (int i=0; i<NCELLS; i++){
		xpin0[NDIM*i] = 0.5*L + 0.6*sqrt(meso2Dobj.geta0(0))*cos(th);
		xpin0[NDIM*i + 1] = 0.5*L + 0.6*sqrt(meso2Dobj.geta0(0))*sin(th);
		th += (2.0*PI)/NCELLS;
	}

	// draw pins to box center
	meso2Dobj.mesoPinFIRE(xpin0, Ftol, dt0, kcspring);

	// set aging parameters
	meso2Dobj.setbetaEff(betaEff);
	meso2Dobj.setcL(cL);
	meso2Dobj.setaL(aL);
	meso2Dobj.setcB(cB);
	meso2Dobj.setcKb(cKb);
	meso2Dobj.setkbi(kb0);

	// initialize adhesive network contacts
	meso2Dobj.initializeMesophyllBondNetwork();

	// run stretching simulation to create network
	meso2Dobj.mesoPinExtension(Ftol, dt0, hmax, dh, dhprint, 2.0*kcspring);


	// say goodbye
	cout << "\n** Finished mesoTriplet2D.cpp, ending. " << endl;

	return 0;
}