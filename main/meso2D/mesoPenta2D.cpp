// FILE to create mini-network with 6 (5 pullers, 1 center) mesophyll cells in 2D
// 
// * Cells are pulled to center to start
// * Contacts are formed
// * As cells are pulled away from center (no size change), shapes rigidify with aging perimeters, angles + kb
// * network is created until hmax is reached
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/meso2D/mesoPenta2D.cpp src/*.cpp -o meso.o
// ./meso.o 24 1.06 1e-3 1e-6 10.0 1.0 0.5 0.01 1.0 0.001 0 1 pos.test
// 
// 
// Parameter input list
// 1. n1: 				number of vertices on first particle
// 2. calA0: 			preferred initial shape parameter for all particles
// 3. dh 				step size
// 4. kb0 				initial bending energy
// 5. betaEff: 			effective temperature, sets contact breaking
// 6. ctcdel: 			amount of contact-dependent adhesion (0=none, 1=all)
// 7. ctch: 			dimensionless bond breaking energy (default=0.5)
// 8. cL: 				perimeter aging parameter
// 9. aL: 				distribution of aging parameter to contact (0) vs void (1)
// 10. cB: 				preferred angle aging parameter
// 11. cKb; 			bending energy aging parameter
// 12. seed: 			seed for random number generator
// 13. positionFile: 	string of path to output file with position/configuration data
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
const int NCELLS 				= 6;		// always 6 cells (5 boundary, 1 center)
const double phi0 				= 0.1;		// initial packing fraction, for viz
const double hmax 				= 2.5;		// max step length
const double dhprint 			= 0.01;		// dh before print step
const double boxLengthScale 	= 2.5;		// neighbor list box size in units of initial l0
const double dt0 				= 1e-2;		// initial magnitude of time step in units of MD time
const double Ftol 				= 1e-12; 	// force tolerance
const double kcspring 			= 1.0; 		// spring connecting to centers
const double kl 				= 1.0; 		// perimeter spring constant

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int n1, NVMAX, seed;
	double calA0, betaEff, cL, aL, cB, cKb, L, dh, kb0, ctcdel, ctch, rtmp;

	// read in parameters from command line input
	string n1_str 			= argv[1];
	string calA0_str 		= argv[2];
	string dh_str 			= argv[3];
	string kb0_str 			= argv[4];
	string betaEff_str 		= argv[5];
	string ctcdel_str 		= argv[6];
	string ctch_str 		= argv[7];
	string cL_str 			= argv[8];
	string aL_str 			= argv[9];
	string cB_str 			= argv[10];
	string cKb_str 			= argv[11];
	string seed_str 		= argv[12];
	string positionFile 	= argv[13];
	string bondFile 		= argv[14];

	// using sstreams to get parameters
	stringstream n1ss(n1_str);
	stringstream calA0ss(calA0_str);
	stringstream dhss(dh_str);
	stringstream kb0ss(kb0_str);
	stringstream betaEffss(betaEff_str);
	stringstream ctcdelss(ctcdel_str);
	stringstream ctchss(ctch_str);
	stringstream cLss(cL_str);
	stringstream aLss(aL_str);
	stringstream cBss(cB_str);
	stringstream cKbss(cKb_str);
	stringstream seedss(seed_str);

	// read into data
	n1ss 			>> n1;
	calA0ss 		>> calA0;
	dhss			>> dh;
	kb0ss 			>> kb0;
	betaEffss 		>> betaEff;
	ctcdelss  		>> ctcdel;
	ctchss  		>> ctch;
	cLss 			>> cL;
	aLss 			>> aL;
	cBss 			>> cB;
	cKbss 			>> cKb;
	seedss 			>> seed;

	// check inputs
	if (ctcdel < 0.0 || ctcdel > 1.0){
		cout << "** ERROR: ctcdel = " << ctcdel << ", but it needs to be between 0 and 1. Ending here. " << endl;
		return 1;
	}

	// instantiate object
	meso2D meso2Dobj(NCELLS, seed);

	// SET PBC -> 0, NO NEED FOR TRIPLET CELLS
	meso2Dobj.setpbc(0,false);
	meso2Dobj.setpbc(1,false);
	meso2Dobj.setkl(kl);

	// open position config file
	meso2Dobj.openPosObject(positionFile);
	meso2Dobj.openCTCObject(bondFile);

	// initialize particles
	meso2Dobj.initializeMesophyllCells(0.0, calA0, phi0, Ftol, n1);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// set max # of vertices
	NVMAX = 2*meso2Dobj.getNVTOT();
	meso2Dobj.setNVMAX(NVMAX);

	// put initial pins in box center
	L = meso2Dobj.getL(0);
	vector<double> xpin0(NDIM*NCELLS,0.0);
	double th = 0.0;
	xpin0[0] = 0.5*L;
	xpin0[1] = 0.5*L;
	for (int i=1; i<NCELLS; i++){
		rtmp = sqrt(meso2Dobj.geta0(0)/PI) + meso2Dobj.getl0(0);
		xpin0[NDIM*i] = 0.5*L + 1.5*rtmp*cos(th);
		xpin0[NDIM*i + 1] = 0.5*L + 1.5*rtmp*sin(th);
		th += (2.0*PI)/(NCELLS-1);
	}

	// set aging parameters
	meso2Dobj.setbetaEff(betaEff);
	meso2Dobj.setctcdel(ctcdel);
	meso2Dobj.setctch(ctch);
	meso2Dobj.setcL(cL);
	meso2Dobj.setaL(aL);
	meso2Dobj.setcB(cB);
	meso2Dobj.setcKb(cKb);
	meso2Dobj.setkbi(kb0);

	// draw pins to box center
	meso2Dobj.setkl(kl);
	meso2Dobj.mesoPinFIRE(xpin0, Ftol, dt0, 0.1*kcspring);

	// initialize adhesive network contacts
	meso2Dobj.initializeMesophyllBondNetwork();

	// run stretching simulation to create network
	meso2Dobj.mesoPinExtension(Ftol, dt0, hmax, dh, dhprint, kcspring, 0);

	// say goodbye
	cout << "\n** Finished mesoPenta2D.cpp, ending. " << endl;

	return 0;
}