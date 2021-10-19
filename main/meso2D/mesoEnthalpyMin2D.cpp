// FILE to create network of mesophyll cells in 2D
//
// * Cells begin at phi0 with no bending energy
// * Compressed to target pressure
// * Cells pulled away WITH MINIMIZED ENTHALPY!!!
// * network is created until phimin is reached
//
// Compilation command:
// g++ -O3 --std=c++11 -I src main/meso2D/mesoEnthalpyMin2D.cpp src/*.cpp -o meso.o
// ./meso.o 16 24 1.12 1e-4 5 1e-3 2 1e-5 0.25 2 0.5 1 pos.test
//
//
// Parameter input list
// 1. NCELLS: 			number of cells
// 2. n1: 				number of vertices on first particle
// 3. calA0: 			preferred initial shape parameter for all particles
// 4. kb0: 				initial amount of bending energy
// 5. betaEff: 			inverse bond breaking temperature
// 6. da0: 				particle growth rate
// 7. dl0: 				perimeter growth rate
// 8. P0: 				fixed pressure
// 9. ctch 				breaking strength
// 10. cL: 				perimeter aging parameter
// 11. cB: 				preferred angle aging parameter
// 12. seed: 			seed for random number generator
// 13. positionFile: 	string of path to output file with position/configuration data

// header files
#include "meso2D.h"
#include <sstream>

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

// global constants
const double plotCompression = 0;  	// whether or not to plot compression
const double dphiGrow = 0.01;	   	// packing fraction increment during initial growth step
const double delShrink = 1e-3;		// fractional change in effective box length during extension
const double boxLengthScale = 3.0; 	// neighbor list box size in units of initial l0
const double phi0 = 0.5;		   	// initial packing fraction
const double dt0 = 1e-2;		   	// initial magnitude of time step in units of MD time
const double Ptol = 1e-6;		   	// target pressure in initial compression
const double Ftol = 1e-10; 			// force tolerance
const double dPtol = 1e-10; 		// fixed pressure tolerance
const double phiMin = 0.4;			// minimum packing fraction in decompression algorithm
const double T0 = 1e-2; 			// temperature for jamming preparation protocol
const double trun = 50.0; 			// amount of time to run annealing
const double kl = 0.1; 				// perimeter spring constant
const double kc = 0.5; 				// interaction spring constant
const double dispersion = 0.1; 		// polydispersity (fixed)
const double aL = 1.0; 				// distribute aging to void only
const double cKb = 0.0; 			// bending rigidity aging
const int NMINSKIP = 10; 			// number of frames to skip

// set parameters
const double ctcdel = 1.0;

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int NCELLS, n1, seed;
	double dispersion, calA0, kb0, betaEff, P0, da0, dl0, ctch, cL, cB, NVMAX;

	// read in parameters from command line input
	string NCELLS_str 		= argv[1];
	string n1_str 			= argv[2];
	string calA0_str 		= argv[3];
	string kb0_str 			= argv[4];
	string betaEff_str 		= argv[5];
	string da0_str 			= argv[6];
	string dl0_str 			= argv[7];
	string P0_str 			= argv[8];
	string ctch_str 		= argv[9];
	string cL_str 			= argv[10];
	string cB_str 			= argv[11];
	string seed_str 		= argv[12];
	string positionFile 	= argv[13];

	// using sstreams to get parameters
	stringstream NCELLSss(NCELLS_str);
	stringstream n1ss(n1_str);
	stringstream calA0ss(calA0_str);
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
	NCELLSss >> NCELLS;
	n1ss >> n1;
	calA0ss >> calA0;
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
	meso2D meso2Dobj(NCELLS, seed);

	// open position config file
	meso2Dobj.openPosObject(positionFile);

	// initialize particles are bidisperse
	meso2Dobj.initializeMesophyllCells(dispersion, calA0, phi0, Ftol, n1);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// set max # of vertices
	NVMAX = 5*meso2Dobj.getNVTOT();
	meso2Dobj.setNVMAX(NVMAX);

	// jam to target pressure
	meso2Dobj.vertexAnneal2Jam2D(&dpm::repulsiveForceUpdate, Ftol, Ptol, dt0, dphiGrow, T0, trun, plotCompression);

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
	meso2Dobj.printMesoNetworkCTCS2D();

	// run stretching simulation to create network
	meso2Dobj.mesoNetworkEnthalpyMin(&meso2D::mesoNetworkForceUpdate, Ftol, dPtol, dt0, da0, dl0, P0, phiMin, NMINSKIP);

	// say goodbye
	cout << "\n** Finished mesoNetwork2D.cpp, ending. " << endl;

	return 0;
}