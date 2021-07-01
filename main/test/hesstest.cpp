// File to compute dynamical matrix of jammed configuration of DPM particles at target pressure
// 
// Will:
// 		create bidisperse DPM particles, 
// 		set constants,
// 		place particle centers, 
// 		relax shapes + positions, 
// 		compress to target pressure Ptol (jamming is limit of Ptol -> 0), 
// 		** COMPUTE HESSIAN, 
// 		print configuration
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/hesstest.cpp src/*.cpp -o test.o
// 
// Run command:
// ./test.o


// header files
#include "dpm.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	bool plotCompression = 0;
	int NCELLS = 12, nsmall = 16, seed = 1;
	double phi0 = 0.6, calA0 = 1.1, smallfrac = 0.5, sizefrac = 1.4, Ftol = 1e-12, Ptol = 1e-8, dt0 = 1e-2;
	double ka = 1.0, kl = 1.0, kb = 0, kc = 1.0, thA = 12.0, thK = 3.0, boxLengthScale = 2.5;

	// name of output file
	string posf = "pos.test";
	string hessf = "hess.test";

	// instantiate object
	dpm configobj2D(NCELLS, NDIM, seed);

	// open position config file
	configobj2D.openPosObject(posf);

	// open hessian ofstream object
	ofstream hessout(hessf.c_str());
	if (!hessout.is_open()){
		cerr << "** Hessian file string " << hessf << " did not open, ending here." << endl;
		return 1;
	}

	// set spring constants
	configobj2D.setka(ka);
	configobj2D.setkl(kl);
	configobj2D.setkb(kb);
	configobj2D.setkc(kc);

	// initialize particles are bidisperse
	configobj2D.bidisperse2D(calA0, nsmall, smallfrac, sizefrac);

	// set preferred angle to sinusoidal
	// configobj2D.sinusoidalPreferredAngle(thA, thK);

	// initialize particle positions
	configobj2D.initializePositions2D(phi0, Ftol);

	// initialize neighbor linked list
	configobj2D.initializeNeighborLinkedList2D(boxLengthScale);

	// compress to target packing fraction
	double phi0Target = 1.0, dphi0 = 0.005;
	configobj2D.vertexJamming2D(&dpm::forceUpdate, Ftol,Ptol,dt0,dphi0,plotCompression);

	// compute Hessian
	int vertDOF = configobj2D.getvertDOF();
	Eigen::MatrixXd M(vertDOF,vertDOF);
	Eigen::MatrixXd H(vertDOF,vertDOF);
	Eigen::MatrixXd S(vertDOF,vertDOF);
	configobj2D.dpmHessian2D(H,S);

	// fill dynamical matrix
	for (int k=0; k<vertDOF; k++){
		for (int l=0; l<vertDOF; l++)
			M(k,l) = H(k,l) - S(k,l);
	}

	// print hessian data
	configobj2D.printHessianEigenvalues2D(hessout,M);

	// close hessian output file
	hessout.close();

	// say goodbye
	cout << "\n\n** Finished hesstest.cpp, ending. " << endl;

	return 0;
}