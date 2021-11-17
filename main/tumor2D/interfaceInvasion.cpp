// File to read in initialized interface configuration, run invasion protocol
// 
// Reads in configuration, sets constants, invades
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/tumor2D/interfaceInvasion.cpp src/*.cpp -o tumor.o
// 
// Example execution:
// ./tumor.o tumor_input.test 1e5 1e3 0.01 0.02 0.1 0.1 0.2 0.5 0.1 0.9 0.01 1e-3 1 pos.test
// 
// 


// header files
#include "tumor2D.h"
#include <sstream>

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

// global constants
const double ka = 1.0;					// area force spring constant (should be unit)
const double kl = 0.2; 				// contractility spring constant
const double kb = 1e-3;					// bending energy
const double kc = 0.05;					// interaction force spring constant
const double gamtt = 0.0; 				// surface tension
const double boxLengthScale = 3.0;		// neighbor list box size in units of initial l0
const double dt0 = 2e-2;				// initial magnitude of time step in units of MD time
const double Ftol = 1e-12;

int main(int argc, char const *argv[])
{
	// forces for simulations
	tumor2DMemFn invasionForceUpdate = nullptr;
	invasionForceUpdate = &tumor2D::stickyTumorInterfaceForceUpdate;

	// local variables to be read in
	int NT, NPRINTSKIP, seed;
	double NTdbl, NPRINTSKIPdbl, l1, l2, v0, Dr0, Ds, kecm, ecmbreak, dDr, dPsi, Drmin;

	// read in parameters from command line input
	string inputFile 		= argv[1];				// input file with initial configuration
	string NT_str 			= argv[2];				// # of time steps
	string NPRINTSKIP_str 	= argv[3];				// # of steps between prints
	string l1_str 			= argv[4];				// attraction strength (must be < l2)
	string l2_str 			= argv[5];				// attraction range (must be > l1)
	string v0_str 			= argv[6];				// tumor cell crawling speed
	string Dr0_str 			= argv[7];				// initial angular diffusion
	string Ds_str 			= argv[8];				// spread of velocity around cell perimeter
	string kecm_str 		= argv[9];				// ecm adhesion strength
	string ecmbreak_str 	= argv[10];				// ecm adhesion range
	string dDr_str 			= argv[11];				// step down when cells in range of adipocytes
	string dPsi_str 		= argv[12];				// relaxation toward x-direction
	string Drmin_str 		= argv[13];				// minimum angular diffusion
	string seed_str 		= argv[14];				// seed for rng
	string positionFile 	= argv[15];				// output file string

	// using sstreams to get parameters
	stringstream NTss(NT_str);
	stringstream NPRINTSKIPss(NPRINTSKIP_str);
	stringstream l1ss(l1_str);
	stringstream l2ss(l2_str);
	stringstream v0ss(v0_str);
	stringstream Dr0ss(Dr0_str);
	stringstream Dsss(Ds_str);
	stringstream kecmss(kecm_str);
	stringstream ecmbreakss(ecmbreak_str);
	stringstream dDrss(dDr_str);
	stringstream dPsiss(dPsi_str);
	stringstream Drminss(Drmin_str);
	stringstream seedss(seed_str);

	// read into data
	NTss 			>> NTdbl;
	NPRINTSKIPss 	>> NPRINTSKIPdbl;
	l1ss 			>> l1;
	l2ss 			>> l2;
	v0ss 			>> v0;
	Dr0ss 			>> Dr0;
	Dsss 			>> Ds;
	kecmss			>> kecm;
	ecmbreakss 		>> ecmbreak;
	dDrss 			>> dDr;
	dPsiss 			>> dPsi;
	Drminss  		>> Drmin;
	seedss 			>> seed;

	// cast step dbls to ints
	NT = (int)NTdbl;
	NPRINTSKIP = (int)NPRINTSKIPdbl;

	// instantiate object
	tumor2D tumor2Dobj(inputFile,seed);

	// open position config file
	tumor2Dobj.openPosObject(positionFile);

	// set spring constants
	tumor2Dobj.setka(ka);
	tumor2Dobj.setkl(kl);
	tumor2Dobj.setkb(kb);
	tumor2Dobj.setkc(kc);
	tumor2Dobj.setgamtt(gamtt);

	// activity parameters
	tumor2Dobj.setv0(v0);
	tumor2Dobj.setDr0(Dr0);
	tumor2Dobj.setDs(Ds);
	tumor2Dobj.setkecm(kecm);
	tumor2Dobj.setecmbreak(ecmbreak);

	// attraction parameters
	tumor2Dobj.setl1(l1);
	tumor2Dobj.setl2(l2);

	// time step in MD time units
	tumor2Dobj.setdt(dt0);

	// initialize neighbor linked list
	tumor2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// run FIRE to relax forces fully
	//tumor2Dobj.tumorFIRE(invasionForceUpdate,Ftol,0.2*dt0);

	// invasion
	cout << "Running invasion protocol..." << endl;
	// tumor2Dobj.invasion(invasionForceUpdate,dDr,dPsi,Drmin,NT,NPRINTSKIP);
	tumor2Dobj.invasionConstP(invasionForceUpdate,dDr,dPsi,Drmin,NT,NPRINTSKIP);

	// say goodbye
	cout << "\n** Finished interfaceInvasion.cpp, ending. " << endl;

	return 0;
}
