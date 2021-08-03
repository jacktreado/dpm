// File to run crawling monolayer simulation
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/tumor2D/monolayerCrawling.cpp src/*.cpp -o tumor.o
// 
// Example execution:
// ./tumor.o 16 24 1.16 0.01 0.1 0.01 5.0 1 pos.test
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
const double kl = 1.0; 					// contractility spring constant
const double kb = 0.0;					// bending energy
const double kc = 1.0;					// interaction force spring constant (should be unit)

const double v0 = 0.05; 				// crawling speed 
const double Ds = 0.1;					// velocity spread around boundary
const double l2 = 0.05; 				// attraction range

const double gamtt = 0.0; 				// surface tension
const double boxLengthScale = 2.0;		// neighbor list box size in units of initial l0
const double dt0 = 2e-2;				// initial magnitude of time step in units of MD time

// 
// 
// TO DO:
// -- Finish parameter read in 
// -- Write bash file
// -- Test locally
// -- push to cluster



int main(int argc, char const *argv[])
{
	// local variables to be read in
	int tN, aNV, seed;
	double tCalA0, l1, Dr, gamtt, tau;

	// read in parameters from command line input
	string tN_str 			= argv[1];				// number of tumor cells
	string tNV_str 			= argv[2];				// tumor cell mean vertex number
	string tCalA0_str 		= argv[3];				// tumor cell shape parameter
	string l1_str 			= argv[4];				// attraction strength (must be < 0.5)
	string Dr_str 			= argv[5];				// initial angular diffusion
	string gamtt_str 		= argv[6]; 				// surface tension
	string tau_str 			= argv[7]; 				// vicsek alignment timescale
	string seed_str 		= argv[10];				// seed for rng
	string positionFile 	= argv[11];				// string for postion output file

	// // using sstreams to get parameters
	// stringstream aNss(aN_str);
	// stringstream aNVss(aNV_str);
	// stringstream tNVss(tNV_str);
	// stringstream aDispss(aDisp_str);
	// stringstream tDispss(tDisp_str);
	// stringstream aCalA0ss(aCalA0_str);
	// stringstream tCalA0ss(tCalA0_str);
	// stringstream areaRatioss(areaRatio_str);
	// stringstream prtss(prt_str);
	// stringstream seedss(seed_str);

	// // read into data
	// aNss 			>> aN;
	// aNVss 			>> aNV;
	// tNVss 			>> tNV;
	// aDispss 		>> aDisp;
	// tDispss 		>> tDisp;
	// aCalA0ss 		>> aCalA0;
	// tCalA0ss 		>> tCalA0;
	// areaRatioss 	>> areaRatio;
	// prtss 			>> prt;
	// seedss 			>> seed;

	// instantiate object
	tumor2D tumor2Dobj(NCELLS, seed);

	// set periodic boundaries
	tumor2Dobj.setpbc(0,1);
	tumor2Dobj.setpbc(1,1);

	// open position config file
	tumor2Dobj.openPosObject(positionFile);

	// set spring constants
	tumor2Dobj.setka(ka);
	tumor2Dobj.setkl(kl);
	tumor2Dobj.setkb(kb);
	tumor2Dobj.setkc(kc);


	// activity parameters
	tumor2Dobj.setv0(v0);
	tumor2Dobj.setDr0(Dr);
	tumor2Dobj.setDs(Ds);
	tumor2Dobj.settau(tau);
	tumor2Dobj.setl1(l1);
	tumor2Dobj.setl2(l2);

	// initialize tumor cells
	tumor2Dobj.gaussian2D(tDisp, tCalA0, tNV);

	// initialize particle positions
	tumor2Dobj.initializeTumorMonolayerPositions(phi0, Ftol, kwell);

	// initialize neighbor linked list
	tumor2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// force updates for crawling
	crawlingForceUpdate = &tumor2D::stickyTumorForceUpdate;
	psiUpdate = &tumor2D::psiVicsek;

	// relax positions using tumor FIRE
	tumor2Dobj.tumorFIRE(crawlingForceUpdate, Ftol, dt0);

	// -- crawling simming
	tumor2Dobj.setgamtt(gamtt);
	tumor2Dobj.setdt(dt0);
	tumor2Dobj.crawling(crawlingForceUpdate,psiUpdate,NT,NPRINTSKIP);

	// say goodbye
	cout << "\n\n** Finished monolayerCrawling.cpp, ending. " << endl;

	return 0;
}