// File to jam configuration of DPM particles to target pressure
// ** Will specifically jam particles with sinusoidal preferred curvature
// 
// Create bidisperse DPM particles, set constants, place particle centers, 
// relax shapes + positions, compress to target pressure Ptol (jamming is limit of Ptol -> 0), print configuration
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/jam/bidisperseSinusoidalParticleJamming.cpp src/*.cpp -o sinejam.o
// 
// 
// Parameter input list
// 1. NCELLS: 			number of particles
// 2. nsmall: 			number of vertices on small particles (larger particles set by 1.4:1.0 bidispersity)
// 3. calA0: 			preferred shape parameter for all particles
// 4. kl: 				perimeter spring constant
// 5. kb: 				bending spring constant
// 6. thA: 				sinusoidal preferred angle amplitude 
// 7. thK: 				sinusoidal preferred angle wavenumber
// 8. Ptol: 			pressure tolerance, sets distance to jamming
// 9. Ftol: 			force tolerance, sets distance to each energy minimum
// 10. seed: 			seed for random number generator
// 11. positionFile: 	string of path to output file with position/configuration data


// header files
#include "dpm.h"
#include <sstream>

// preprocessor macros
#define NDIM 2
	
// namspace
using namespace std;

// global constants
const bool plotCompression = 0;			// whether or not to plot configuration during jamming (0 saves memory)
const double dphi0 = 0.005;				// packing fraction increment
const double ka = 1.0;					// area force spring constant (should be unit)
const double kc = 1.0;					// interaction force spring constant (should be unit)
const double boxLengthScale = 2.5;		// neighbor list box size in units of initial l0
const double phi0 = 0.5;				// initial packing fraction
const double smallfrac = 0.5;			// fraction of small particles
const double sizeratio = 1.4;			// size ratio between small and large particles
const double dt0 = 1e-2;				// initial magnitude of time step in units of MD time

int main(int argc, char const *argv[])
{
	// local variables to be read in
	int NCELLS, nsmall, seed;
	double calA0, kl, kb, thA, thK, Ptol, Ftol;

	// read in parameters from command line input
	string NCELLS_str 		= argv[1];
	string nsmall_str 		= argv[2];
	string calA0_str 		= argv[3];
	string kl_str 			= argv[4];
	string kb_str 			= argv[5];
	string thA_str 			= argv[6];
	string thK_str 			= argv[7];
	string Ptol_str 		= argv[8];
	string Ftol_str 		= argv[9];
	string seed_str 		= argv[10];
	string positionFile 	= argv[11];

	// using sstreams to get parameters
	stringstream NCELLSss(NCELLS_str);
	stringstream nsmallss(nsmall_str);
	stringstream calA0ss(calA0_str);
	stringstream klss(kl_str);
	stringstream kbss(kb_str);
	stringstream thAss(thA_str);
	stringstream thKss(thK_str);
	stringstream Ptolss(Ptol_str);
	stringstream Ftolss(Ftol_str);
	stringstream seedss(seed_str);

	// read into data
	NCELLSss 		>> NCELLS;
	nsmallss 		>> nsmall;
	calA0ss 		>> calA0;
	klss 			>> kl;
	kbss 			>> kb;
	thAss 			>> thA;
	thKss 			>> thK;
	Ptolss			>> Ptol;
	Ftolss 			>> Ftol;
	seedss 			>> seed;

	// instantiate object
	dpm configobj2D(NCELLS, NDIM, seed);

	// open position config file
	configobj2D.openPosObject(positionFile);

	// set spring constants
	configobj2D.setka(ka);
	configobj2D.setkl(kl);
	configobj2D.setkb(kb);
	configobj2D.setkc(kc);

	// initialize particles are bidisperse
	configobj2D.bidisperse2D(calA0, nsmall, smallfrac, sizeratio);

	// set preferred angle to sinusoidal
	configobj2D.sinusoidalPreferredAngle(thA, thK);

	// initialize particle positions
	configobj2D.initializePositions2D(phi0, Ftol);

	// initialize neighbor linked list
	configobj2D.initializeNeighborLinkedList2D(boxLengthScale);

	// compress to target packing fraction
	configobj2D.vertexJamming2D(Ftol,Ptol,dt0,dphi0,plotCompression);

	// say goodbye
	cout << "\n** Finished bidisperseSinusoidalParticleJamming.cpp, ending. " << endl;

	return 0;
}