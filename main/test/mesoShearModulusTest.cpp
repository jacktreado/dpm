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

	// compute shear modulus numerically, using 1st and 2nd derivatives
	int NVTOT = meso2Dobj.getNVTOT();
	int NVVCTS = 0.5*NVTOT*(NVTOT-1);
	vector<bool> gijtmp(NVVCTS,0);

	// initialize temporary contact network
	meso2Dobj.getMesoVVContactNetwork(gijtmp);






	// -- Compute shear modulus numerically

	// shear strain
	double dgamma = 2e-8;
	double gamma = 0.0;
	int NGAMMA = 50;
	int k = 0;

	// save shear stress and potential energy
	vector<double> UList(NGAMMA+1,0.0);
	vector<double> sxyList(NGAMMA+1,0.0);

	// save initial shear stress and potential energy
	sxyList.at(0) = meso2Dobj.getstress(2);
	UList.at(0) = meso2Dobj.getU();

	// loop over shear strains gamma, relax using FIRE, compute change in Sxy
	for (k=0; k<NGAMMA; k++){
		// update gamma for this iteration
		gamma += dgamma;

		// relax at fixed shear strain + volume, FIXED CONTACT NETWORK
		meso2Dobj.mesoShearStrainFIRE(gamma, Ftol, dt0, gijtmp);

		// save shear stress
		sxyList.at(k+1) = meso2Dobj.getstress(2);

		// save total potential energy
		UList.at(k+1) = meso2Dobj.getU();

		// print
		meso2Dobj.printMesoShearConfigCTCS2D(gamma);
	}





	// compute shear moduli using second derivatives
	vector<double> G_sxy_list(NGAMMA-1,0.0);
	vector<double> G_U_list(NGAMMA-1,0.0);
	vector<double> sxy_U_list(NGAMMA-1,0.0);
	for (k=0; k<NGAMMA-1; k++){
		G_sxy_list.at(k) = -0.5*(sxyList.at(k+2) - sxyList.at(k))/dgamma;
		G_U_list.at(k) = (UList.at(k+2) + UList.at(k) - 2.0*UList.at(k+1))/(dgamma*dgamma);
		sxy_U_list.at(k) = -0.5*(UList.at(k+2) - UList.at(k))/dgamma;
	}


	// print
	cout << "shear modulus G computed by 2 methods:" << endl;
	gamma = dgamma;
	for (k=0; k<NGAMMA-1; k++){
		cout << "k = " << k << ", gamma = " << gamma << ";   sxy = " << sxyList.at(k) << ",  sxy_U = " << sxy_U_list.at(k) << ";   G_sxy = " << G_sxy_list.at(k) << ",  G_U = " << G_U_list.at(k) << endl;
		gamma += dgamma;
	}


	// say goodbye
	cout << "\n** Finished mesoShearModulusTest.cpp, ending. " << endl;

	return 0;
}