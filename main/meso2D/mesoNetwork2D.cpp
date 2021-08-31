// FILE to create network of mesophyll cells in 2D
//
// * Cells begin at phi0 with no bending energy
// * Compressed to target pressure
// * As cells are pulled away, shapes rigidify with aging perimeters, angles + kb
// * network is created until phimin is reached
//
// Compilation command:
// g++ -O3 --std=c++11 -I src main/meso2D/mesoNetwork2D.cpp src/*.cpp -o meso.o
// ./meso.o 16 24 0.1 1.06 10.0 0.01 1 1e-2 1e-6 1 pos.test
//
//
// Parameter input list
// 1. NCELLS: 			number of cells
// 2. n1: 				number of vertices on first particle
// 3. dispersion: 		polydispersity
// 4. calA0: 			preferred initial shape parameter for all particles
// 5. betaEff: 			inverse bond breaking temperature
// 6. cL: 				perimeter aging parameter
// 7. aL: 				distribution of aging to either contacts (0) or void (1)
// 8. cB: 				preferred angle aging parameter
// 9. cKb; 				bending energy aging parameter
// 10. seed: 			seed for random number generator
// 11. positionFile: 	string of path to output file with position/configuration data

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
const double delShrink = 0.001;	   	// fractional change in effective box length during extension
const double dphiPrint = 0.01;	   	// packing fractions to skip between print steps
const double boxLengthScale = 2.5; 	// neighbor list box size in units of initial l0
const double phi0 = 0.5;		   	// initial packing fraction
const double dt0 = 1e-2;		   	// initial magnitude of time step in units of MD time
const double Ptol = 1e-6;		   	// target pressure in initial compression
const double Ftol = 1e-10; 			// force tolerance
const double phiMin = 0.4;			// minimum packing fraction in decompression algorithm
const double T0 = 1e-3; 			// temperature for jamming preparation protocol
const double trun = 50.0; 			// amount of time to run

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int NCELLS, n1, seed;
	double dispersion, calA0, phiMin, betaEff, cL, aL, cB, cKb;

	// read in parameters from command line input
	string NCELLS_str 		= argv[1];
	string n1_str 			= argv[2];
	string disp_str 		= argv[3];
	string calA0_str 		= argv[4];
	string betaEff_str 		= argv[5];
	string cL_str 			= argv[6];
	string aL_str 			= argv[7];
	string cB_str 			= argv[8];
	string cKb_str 			= argv[9];
	string seed_str 		= argv[10];
	string positionFile 	= argv[11];

	// using sstreams to get parameters
	stringstream NCELLSss(NCELLS_str);
	stringstream n1ss(n1_str);
	stringstream dispss(disp_str);
	stringstream calA0ss(calA0_str);
	stringstream betaEffss(betaEff_str);
	stringstream cLss(cL_str);
	stringstream aLss(aL_str);
	stringstream cBss(cB_str);
	stringstream cKbss(cKb_str);
	stringstream seedss(seed_str);

	// read into data
	NCELLSss >> NCELLS;
	n1ss >> n1;
	dispss >> dispersion;
	calA0ss >> calA0;
	betaEffss >> betaEff;
	cLss >> cL;
	aLss >> aL;
	cBss >> cB;
	cKbss >> cKb;
	seedss >> seed;

	// instantiate object
	meso2D meso2Dobj(NCELLS, seed);

	// open position config file
	meso2Dobj.openPosObject(positionFile);

	// initialize particles are bidisperse
	meso2Dobj.initializeMesophyllCells(dispersion, calA0, phi0, Ftol, n1);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// compress to target packing pressure
	double Ftoltmp = Ftol;
	if (Ftol < 1e-10)
		Ftoltmp = Ftol * 100;

	// jam to target pressure
	// meso2Dobj.vertexJamming2D(&dpm::repulsiveForceUpdate, Ftoltmp, Ptol, dt0, dphiGrow, plotCompression);
	meso2Dobj.vertexAnneal2Jam2D(&dpm::repulsiveForceUpdate, Ftoltmp, Ptol, dt0, dphiGrow, T0, trun, plotCompression);

	// set aging parameters
	meso2Dobj.setbetaEff(betaEff);
	meso2Dobj.setcL(cL);
	meso2Dobj.setaL(aL);
	meso2Dobj.setcB(cB);
	meso2Dobj.setcKb(cKb);

	// relax configuration using network + bending
	meso2Dobj.initializeMesophyllBondNetwork();
	meso2Dobj.t0ToCurrent();
	meso2Dobj.mesoFIRE(&meso2D::mesoNetworkForceUpdate, Ftol, dt0);
	meso2Dobj.printMesoNetwork2D();

	// run stretching simulation to create network
	meso2Dobj.mesoNetworkExtension(&meso2D::mesoNetworkForceUpdate, Ftol, dt0, delShrink, dphiPrint, phiMin);

	// say goodbye
	cout << "\n** Finished mesoNetwork2D.cpp, ending. " << endl;

	return 0;
}