// FILE to test void bubble code

// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/mesoBubbleTest.cpp src/*.cpp -o meso.o
// ./meso.o meso_n16.bbin pos.test
//
//
// Parameter input list
// 1. inputFile			file to 
// 2. positionFile: 	string of path to output file with position/configuration data

// header files
#include "meso2D.h"
#include <sstream>

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

// global constants
const double delShrink = 2e-3;		// fractional change in effective box length during extension
const double dphiPrint = 0.01;	   	// packing fractions to skip between print steps
const double boxLengthScale = 2.5; 	// neighbor list box size in units of initial l0
const double phi0 = 0.5;		   	// initial packing fraction
const double dt0 = 1e-2;		   	// initial magnitude of time step in units of MD time
const double Ftol = 1e-10; 			// force tolerance
const double dPtol = 1e-10;			// pressure change tolerance
const double phiMin = 0.4;			// minimum packing fraction in decompression algorithm
const double kl = 0.5; 				// perimeter spring constant
const double kc = 0.5; 				// interaction spring constant

// set parameters
const double ctch = 0.5;
const double cKb = 0;


int main(int argc, char const *argv[])
{
	// local variables to be read in
	int NVMAX;;

	// read in parameters from command line input
	string inputFile 		= argv[1];
	string positionFile 	= argv[2];

	// instantiate object
	meso2D meso2Dobj(inputFile,1.0,1);

	// open position config file
	meso2Dobj.openPosObject(positionFile);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// set max # of vertices
	NVMAX = 5*meso2Dobj.getNVTOT();
	meso2Dobj.setNVMAX(NVMAX);

	// set aging parameters
	meso2Dobj.setbetaEff(50.0);
	meso2Dobj.setctcdel(1.0);
	meso2Dobj.setctch(0.5);
	meso2Dobj.setcL(0.0);
	meso2Dobj.setaL(0.0);
	meso2Dobj.setcB(0.0);
	meso2Dobj.setcKb(0.0);
	meso2Dobj.setkbi(0.01);
	meso2Dobj.setkl(1.0);
	meso2Dobj.setkc(1.0);

	// set bonds, preferred curvatures
	meso2Dobj.initializeMesoBubbleBondNetwork();
	meso2Dobj.t0ToCurrent();
	// meso2Dobj.printMesoNetwork2D();

	// run void bubble growth sim
	meso2Dobj.mesoBubbleEnthalpyMin(&meso2D::mesoNetworkForceUpdate, Ftol, dPtol, dt0, 0.01, 1e-4, 0.3, 1);

	// say goodbye
	cout << "\n** Finished mesoBubbleTest.cpp, ending. " << endl;

	return 0;
}