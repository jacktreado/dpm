// header files
#include "tumor2D.h"
#include <sstream>

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

// global constants
const double dphi0= 0.01;	   		// packing fraction increment during initial growth step
const double boxLengthScale = 2.0; 	// neighbor list box size in units of initial l0
const double phi0 = 0.5;		   	// initial packing fraction
const double dt0 = 1e-2;		   	// initial magnitude of time step in units of MD time
const double Ftol = 1e-10;			// force tolerance during energy min
const double Ptol = 1e-6;		   	// target pressure in initial compression


// mechanical constants
const double ka = 1.0;
const double kl = 1.0;
const double kb = 0.0;
const double kc = 1.0;

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int NCELLS, aN, tN, aNV, tNV, seed;
	double aDisp, tDisp, aCalA0, tCalA0, areaRatio, prt;

	// read in parameters from command line input
	string aN_str 			= argv[1];
	string aNV_str 			= argv[2];
	string tNV_str 			= argv[3];
	string aDisp_str 		= argv[4];
	string tDisp_str 		= argv[5];
	string aCalA0_str 		= argv[6];
	string tCalA0_str 		= argv[7];
	string areaRatio_str 	= argv[8];
	string prt_str 			= argv[9];
	string seed_str 		= argv[10];
	string positionFile 	= argv[11];

	// using sstreams to get parameters
	stringstream aNss(aN_str);
	stringstream aNVss(aNV_str);
	stringstream tNVss(tNV_str);
	stringstream aDispss(aDisp_str);
	stringstream tDispss(tDisp_str);
	stringstream aCalA0ss(aCalA0_str);
	stringstream tCalA0ss(tCalA0_str);
	stringstream areaRatioss(areaRatio_str);
	stringstream prtss(prt_str);
	stringstream seedss(seed_str);

	// read into data
	aNss 			>> aN;
	aNVss 			>> aNV;
	tNVss 			>> tNV;
	aDispss 		>> aDisp;
	tDispss 		>> tDisp;
	aCalA0ss 		>> aCalA0;
	tCalA0ss 		>> tCalA0;
	areaRatioss 	>> areaRatio;
	prtss 			>> prt;
	seedss 			>> seed;

	// determine number of tumor cells based on areaRatio and prt
	tN = round(aN * areaRatio * (prt/(1.0 - prt)));
	NCELLS = tN + aN;

	// instantiate object
	tumor2D tumor2Dobj(NCELLS, tN, seed);

	// open position config file
	tumor2Dobj.openPosObject(positionFile);

	// set spring constants
	tumor2Dobj.setka(ka);
	tumor2Dobj.setkl(kl);
	tumor2Dobj.setkb(kb);
	tumor2Dobj.setkc(kc);

	// initialize adipocyte and tumor cells
	tumor2Dobj.initializeTumorInterface(aCalA0, tCalA0, aDisp, tDisp, areaRatio, aNV, tNV);

	// initialize particle positions
	tumor2Dobj.initializeTumorInterfacePositions(phi0, Ftol, prt);

	// initialize neighbor linked list
	tumor2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// -- compression to initial condition
	tumor2Dobj.tumorCompression(Ftol,Ptol,dt0,dphi0);

	// print interface
	tumor2Dobj.printTumorInterface(0.0);


	// end
	cout << "interfaceCompression.cpp completed, interface printed to " << positionFile << ", ending. " << endl;
	return 0;
}

