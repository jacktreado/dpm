// FILE to create mini-network with 6 cells in free boundaries
// compilation:
// g++ --std=c++11 -O3 -I src main/test/mesoFreeGrowthTest.cpp src/*.cpp -o test.o

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
const double hmax 				= 1.25;		// max step length
const double dhprint 			= 0.01;		// dh before print step
const double boxLengthScale 	= 4.0;		// neighbor list box size in units of initial l0
const double dt0 				= 2e-2;		// initial magnitude of time step in units of MD time
const double Ftol 				= 1e-12; 	// force tolerance
const double kcspring 			= 1.0; 		// spring connecting to centers
const double kl 				= 0.1; 		// perimeter spring constant
const double kc 				= 0.1; 		// interaction spring constant

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int n1, NVMAX, seed;
	double calA0, betaEff, cL, aL, cB, cKb, L, kb0, ctcdel, ctch, dl0, dphi0, rtmp, a0max, dphiPrint;

	string calA0_str 		= argv[1];
	string kb0_str 			= argv[2];
	string betaEff_str 		= argv[3];
	string cL_str 			= argv[4];
	string dphi0_str 		= argv[5];
	string dl0_str 			= argv[6];

	stringstream calA0ss(calA0_str);
	stringstream kb0ss(kb0_str);
	stringstream betaEffss(betaEff_str);
	stringstream cLss(cL_str);
	stringstream dphi0ss(dphi0_str);
	stringstream dl0ss(dl0_str);

	calA0ss >> calA0;
	kb0ss >> kb0;
	betaEffss >> betaEff;
	cLss >> cL;
	dphi0ss >> dphi0;
	dl0ss >> dl0;

	// parameters
	n1 = 32;
	ctcdel = 1.0;
	ctch = 0.5;
	aL = 1.0;
	cB = 1.0;
	cKb = 0.0;
	dphiPrint = 0.002;
	a0max = 4.0;
	seed = 1;


	// output strings
	string positionFile = "pos.test";
	string bondFile = "bond.test";

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
	meso2Dobj.setkc(kc);

	// open position config file
	meso2Dobj.openPosObject(positionFile);
	meso2Dobj.openCTCObject(bondFile);

	// initialize particles
	meso2Dobj.initializeMesophyllCells(0.0, calA0, phi0, Ftol, n1);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// set max # of vertices
	NVMAX = 5*meso2Dobj.getNVTOT();
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

	// draw pins to box center
	meso2Dobj.setkl(kl);
	meso2Dobj.mesoPinFIRE(xpin0, Ftol, dt0, 0.1*kcspring);
	meso2Dobj.setkbi(kb0);
	meso2Dobj.t0ToCurrent();

	// initialize adhesive network contacts
	meso2Dobj.initializeMesophyllBondNetwork();

	// run stretching simulation to create network
	meso2Dobj.mesoFreeGrowth(&meso2D::mesoNetworkForceUpdate, Ftol, dt0, dl0, dphi0, dphiPrint, a0max);

	// say goodbye
	cout << "\n** Finished mesoFreeGrowthTest.cpp, ending. " << endl;

	return 0;
}