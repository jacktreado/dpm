// FILE to generate initial condition of packing of mesophyll cells in 2D
//
// Compilation command:
// g++ -O3 --std=c++11 -I src main/meso2D/generateMesoInput2D.cpp src/*.cpp -o meso.o
// ./meso.o 16 24 1.12 1 meso.input
//
//
// Parameter input list
// 1. NCELLS: 			number of cells
// 2. n1: 				number of vertices on first particle
// 3. calA0: 			preferred initial shape parameter for all particles
// 5. seed: 			seed for random number generator
// 6. positionFile: 	string of path to output file with position/configuration data

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
const double boxLengthScale = 3.0; 	// neighbor list box size in units of initial l0
const double phi0 = 0.5;		   	// initial packing fraction
const double dt0 = 1e-2;		   	// initial magnitude of time step in units of MD time
const double Ptol = 1e-7;		   	// target pressure in initial compression
const double Ftol = 1e-12; 			// force tolerance
const double T0 = 1e-2; 			// temperature for jamming preparation protocol
const double trun = 50.0; 			// amount of time to run annealing
const double kl = 0.5; 				// perimeter spring constant
const double kc = 0.5; 				// interaction spring constant
const double dispersion = 0.1; 		// polydispersity (fixed)

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int NCELLS, n1, seed;
	double calA0;

	// read in parameters from command line input
	string NCELLS_str 		= argv[1];
	string n1_str 			= argv[2];
	string calA0_str 		= argv[3];
	string seed_str 		= argv[4];
	string positionFile 	= argv[5];

	// using sstreams to get parameters
	stringstream NCELLSss(NCELLS_str);
	stringstream n1ss(n1_str);
	stringstream calA0ss(calA0_str);
	stringstream seedss(seed_str);

	// read into data
	NCELLSss >> NCELLS;
	n1ss >> n1;
	calA0ss >> calA0;
	seedss >> seed;

	// instantiate object
	meso2D meso2Dobj(NCELLS, seed);

	// open position config file
	meso2Dobj.openPosObject(positionFile);

	// initialize particles are bidisperse
	meso2Dobj.initializeMesophyllCells(dispersion, calA0, phi0, Ftol, n1);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// jam to target pressure
	meso2Dobj.vertexAnneal2Jam2D(&dpm::repulsiveForceUpdate, Ftol, Ptol, dt0, dphiGrow, T0, trun, plotCompression);
	meso2Dobj.printMesoNetwork2D();

	// say goodbye
	cout << endl;
	cout << "** Finished generateMesoInput2D.cpp, saved input to " << positionFile << ". " << endl;
	cout << "** Ending. " << endl;

	return 0;
}