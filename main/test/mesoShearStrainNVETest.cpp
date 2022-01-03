// test main file for meso shear modulus

// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/mesoShearModulusTest.cpp src/*.cpp -o test.o

// header files
#include "meso2D.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int seed = 1;
	double Ftol = 1e-12, Ptol = 1e-6, dt0 = 0.01, dphi0 = 0.01;
	double boxLengthScale = 2.5, betaEff = 50.0, ctcdel = 1.0, ctch = 0.5, cL = 0, aL = 1.0, cB = 0.0, cKb = 0.0, kl = 1.0, kc = 1.0, kb0 = 1e-3;

	// pointer to dpm member function (should pt to null)
	dpmMemFn jammingForceUpdate = nullptr;
	meso2DMemFn mesoForceUpdate = nullptr;

	string inputf = "meso.input";
	string posf = "pos.test";
	string shearf = "shear.test";

	// instantiate object
	meso2D meso2Dobj(inputf,seed);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// open config, hessian and shear file
	meso2Dobj.openPosObject(posf);
	meso2Dobj.openHessObject(shearf);

	// set parameters
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
	meso2Dobj.t0ToCurrent();
	meso2Dobj.printMesoShearConfigCTCS2D(0.0);

	// run NVE protocol at fixed shear strain
	


	// say goodbye
	cout << "\n** Finished mesoShearStrainNVE.cpp, ending. " << endl;

	return 0;
}