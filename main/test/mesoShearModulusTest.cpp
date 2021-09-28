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
	bool plotCompression = 0;
	int NCELLS = 16, n1 = 24, seed = 1;
	double phi0 = 0.7, calA0 = 1.08, dispersion = 0.1, Ftol = 1e-12, Ptol = 1e-5, dt0 = 0.02, dphi0 = 0.01;
	double boxLengthScale = 3.0, betaEff = 1.0, ctcdel = 1.0, ctch = 0.5, cL = 0.01, aL = 1.0, cB = 0.0, cKb = 0.0, kb0 = 1e-3;

	// pointer to dpm member function (should pt to null)
	dpmMemFn jammingForceUpdate = nullptr;
	meso2DMemFn mesoForceUpdate = nullptr;

	string inputf = "meso_input.test";
	string posf = "pos.test";
	string hessf = "hess.test";

	// instantiate object
	meso2D mesoHessObj(inputf,seed);

	// open position config file
	mesoHessObj.openPosObject(posf);

	// initialize neighbor linked list
	mesoHessObj.initializeNeighborLinkedList2D(boxLengthScale);

	// set aging parameters
	mesoHessObj.setbetaEff(betaEff);
	mesoHessObj.setctcdel(ctcdel);
	mesoHessObj.setctch(ctch);
	mesoHessObj.setcL(cL);
	mesoHessObj.setaL(aL);
	mesoHessObj.setcB(cB);
	mesoHessObj.setcKb(cKb);
	mesoHessObj.setkbi(kb0);

	// relax configuration using network + bending
	mesoForceUpdate = &meso2D::mesoNetworkForceUpdate;
	mesoHessObj.initializeMesophyllBondNetwork();
	mesoHessObj.t0ToCurrent();
	mesoHessObj.mesoFIRE(mesoForceUpdate, Ftol, dt0);
	mesoHessObj.printMesoNetwork2D();

	// compute shear modulus using dynamical matrix	
	cout << "** Computing shear modulus using dynamical matrix ... " << endl;
	mesoHessObj.openHessObject(hessf);
	double GDM = mesoHessObj.mesoPrintLinearResponse();

	// get initial shear stress
	cout << "** Now getting shear stress, checking against Lees-Edwards boundary conditions" << endl;
	double sxycurr, sxyold;
	double L = mesoHessObj.getL(0);
	sxyold = mesoHessObj.getstress(2)*(L*L);







	// duplicate system
	cout << "** Duplicating jammed system, computing shear modulus using Lees-Edwards" << endl;
	meso2D mesoShearObj(NCELLS,seed);
	mesoShearObj = mesoHessObj;

	// compute shear modulus via Lees-Edwards boundary conditions
	double dgamma = 1e-8;
	double gamma = dgamma;
	int NGAMMA = 10;
	vector<double> GLE(NGAMMA,0.0);
	for (int g=0; g<NGAMMA; g++){
		// relax using FIRE
		cout << "** In main, minimizing at g = " << g << ", gamma = " << gamma << "...";
		mesoShearObj.mesoShearStrainFIRE(gamma,Ftol,dt0);

		// get out shear stress
		sxycurr = mesoShearObj.getstress(2)*(L*L);
		cout << "sxy = " << sxycurr << endl;

		// compute shear modulus
		GLE.at(g) = -(sxycurr - sxyold)/dgamma;

		// update gamma
		gamma += dgamma;
		sxyold = sxycurr;
	}

	// print results
	cout << "G (dynamical matrix) = " << GDM << endl;
	cout << "G (Lees Edwards, dgamma = " << dgamma << "): " << endl;
	gamma = dgamma;
	for (int g=0; g<NGAMMA; g++){
		cout << "gamma=" << gamma << ";  G=" << GLE.at(g) << endl;
		gamma += dgamma;
	}

	// say goodbye
	cout << "\n** Finished mesoShearModulusTest.cpp, ending. " << endl;

	return 0;
}