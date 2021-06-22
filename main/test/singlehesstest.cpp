// File to compute dynamical matrix of single DPM particle at force balance
// 
// Will:
// 		initialize particle, 
// 		set constants,
// 		relax shapes + positions, 
// 		** COMPUTE HESSIAN, 
// 		print configuration + eigenvalues
// 
// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/singlehesstest.cpp src/*.cpp -o test.o
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
	int NCELLS = 1, n = 16, seed = 1;
	double phi0 = 0.2, calA0 = 1.2, Ftol = 1e-14, Ptol = 1e-8, dt0 = 1e-2;
	double ka = 1.0, kl = 1.0, kb = 0.0, kc = 0.0, thA = 12.0, thK = 3.0, boxLengthScale = 2.5;

	// name of output file
	string posf = "pos.test";
	string hessf = "hess.test";

	// instantiate object
	dpm configobj2D(NCELLS, NDIM, seed);

	// open position config file
	configobj2D.openPosObject(posf);
	configobj2D.openHessObject(hessf);

	// set spring constants
	configobj2D.setka(ka);
	configobj2D.setkl(kl);
	configobj2D.setkb(kb);
	configobj2D.setkc(kc);

	// initialize particles are bidisperse
	configobj2D.gaussian2D(0.0, calA0, n);

	// set preferred angle to sinusoidal
	// configobj2D.sinusoidalPreferredAngle(thA, thK);

	// initialize particle positions
	configobj2D.initializePositions2D(phi0, Ftol);

	// initialize neighbor linked list
	configobj2D.initializeNeighborLinkedList2D(boxLengthScale);

	// use FIRE to relax vertex positions
	configobj2D.vertexFIRE2D(Ftol, dt0);

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

	// print data to file
	configobj2D.printConfiguration2D();
	configobj2D.printMatrixEigenvalues2D(M);


	// say goodbye
	cout << "\n\n** Finished singlehesstest.cpp, ending. " << endl;

	return 0;
}