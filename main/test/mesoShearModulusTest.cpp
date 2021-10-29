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
	int seed = 1;
	double dispersion = 0.1, Ftol = 1e-12, Ptol = 1e-6, dt0 = 0.01, dphi0 = 0.01;
	double boxLengthScale = 10, betaEff = 1.0, ctcdel = 1.0, ctch = 0.5, cL = 0.01, aL = 1.0, cB = 0.0, cKb = 0.0, kb0 = 0;

	// pointer to dpm member function (should pt to null)
	dpmMemFn jammingForceUpdate = nullptr;
	meso2DMemFn mesoForceUpdate = nullptr;

	string inputf = "meso.input";
	string posf = "pos.test";
	string hessf = "hess.test";
	string shearf = "shear.test";

	// instantiate object
	meso2D mesoHessObj(inputf,seed);

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
	// mesoHessObj.initializeMesophyllBondNetwork();
	mesoHessObj.t0ToCurrent();
	mesoHessObj.mesoFIRE(mesoForceUpdate, Ftol, dt0);

	// duplicate system
	cout << "** Duplicating jammed system, computing shear modulus using Lees-Edwards" << endl;
	int NCELLS = mesoHessObj.getNCELLS();
	meso2D mesoShearObj(NCELLS,seed);
	mesoShearObj = mesoHessObj;

	// open position config file to normal system, print for first time
	mesoHessObj.openPosObject(posf);
	mesoHessObj.printMesoNetwork2D();

	// compute shear modulus using dynamical matrix	
	cout << "** Computing shear modulus using dynamical matrix ... " << endl;
	mesoHessObj.openHessObject(hessf);
	double GDM = mesoHessObj.mesoPrintLinearResponse();

	// get initial shear stress
	cout << "** Now getting shear stress, checking against Lees-Edwards boundary conditions" << endl;
	double sxycurr, sxyold, Ucurr, Uold;
	double L = mesoHessObj.getL(0);

	// open print file for sheared system
	mesoShearObj.openPosObject(shearf);
	Uold = mesoHessObj.getU();				// get from Hess obj, should be same as shear obj
	sxyold = mesoShearObj.dUdgamma();

	// compute shear modulus via Lees-Edwards boundary conditions
	double dgamma = 1e-6;
	double gamma = dgamma;
	int NGAMMA = 20;
	vector<double> GLE(NGAMMA,0.0);
	for (int g=0; g<NGAMMA; g++){
		// relax using FIRE
		cout << "** In main, minimizing at g = " << g << ", gamma = " << gamma << "...";
		mesoShearObj.mesoShearStrainFIRE(gamma,Ftol,dt0);
		mesoShearObj.printMesoShearConfig2D(gamma);

		// get out shear stress
		// sxycurr = mesoShearObj.getstress(2)*(L*L);
		Ucurr = mesoShearObj.getU();
		sxycurr = mesoShearObj.dUdgamma();
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