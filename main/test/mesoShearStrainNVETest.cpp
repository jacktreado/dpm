// test main file for meso shear modulus

// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/mesoShearStrainNVETest.cpp src/*.cpp -o test.o
// 
// NOTE: write matlab to read in shear at fixed strain, to visualize NVE with sheared box

// header files
#include "meso2D.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int seed = 1;
	double Ftol = 1e-12, dt0 = 0.02, dphi0 = 0.01;
	double boxLengthScale = 2.5, betaEff = 50.0, ctcdel = 1.0, ctch = 0.5, cL = 0, aL = 1.0, cB = 0.0, cKb = 0.0, kl = 1.0, kc = 1.0, kb0 = 1e-3;

	// pointer to dpm member function (should pt to null)
	dpmMemFn jammingForceUpdate = nullptr;
	meso2DMemFn mesoForceUpdate = nullptr;

	string inputf = "meso.input";
	string posf = "pos.test";
	string enf = "energy.test";

	// instantiate
	meso2D meso2Dobj(inputf,seed);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// open config, hessian and shear file
	meso2Dobj.openPosObject(posf);
	ofstream enout(enf.c_str());
	if (!enout.is_open()){
		cout << "** ERROR: in mesoShearStrainNVETest.cpp, en out file " << enf << " does not open. Ending. " << endl;
		return 1;
	}

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

	// initial vertex-vertex contact network
	int NVTOT = meso2Dobj.getNVTOT();
	int NVVCTS = 0.5*NVTOT*(NVTOT-1);
	vector<bool> gijtmp(NVVCTS,0);

	// initialize contact network
	meso2Dobj.getMesoVVContactNetwork(gijtmp);

	// set gamma and temperature
	double gamma = 0.01;
	double dgamma = 0.001;
	int NGAMMA = 2;
	double T = 5e-2;

	// relax at that gamma
	for (int gg = 0; gg < NGAMMA; gg++){
		meso2Dobj.mesoShearStrainFIRE(gamma, Ftol, dt0, gijtmp);
		gamma += dgamma;
	}

	// run NVE protocol at fixed shear strain
	int NT = 1e6;
	int NPRINTSKIP = 5e3;
	meso2Dobj.mesoShearStrainNVE(enout, gamma, T, dt0, NT, NPRINTSKIP, gijtmp);


	// say goodbye
	cout << "\n** Finished mesoShearStrainNVETest.cpp, ending. " << endl;
	enout.close();
	return 0;
}