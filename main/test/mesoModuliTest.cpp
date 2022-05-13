// test main file for meso shear modulus

// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/mesoModuliTest.cpp src/*.cpp -o test.o

// header files
#include "meso2D.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int seed = 1;
	double Ftol = 1e-12, Ptol = 1e-6, dt0 = 0.01, dphi0 = 0.01;
	double boxLengthScale = 2.5, betaEff = 50.0, ctcdel = 1.0, ctch = 0.5, cL = 0, aL = 1.0, cB = 0.0, cKb = 0.0, kl = 1.0, kc = 1.0, kb0 = 0.0;
	double P0 = 1e-6, dPtol = 1e-10;

	// pointer to dpm member function (should pt to null)
	dpmMemFn jammingForceUpdate = nullptr;
	meso2DMemFn mesoForceUpdate = nullptr;

	string inputf = "meso_n16.input";
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

	// relax configuration at constant pressure using network + bending
	// meso2Dobj.initializeMesophyllBondNetwork();
	meso2Dobj.t0ToCurrent();
	meso2Dobj.mesoFIRE(&meso2D::mesoNetworkForceUpdate, Ftol, dt0);
	meso2Dobj.t0ToCurrent();
	meso2Dobj.mesoEnthalpyFIRE(&meso2D::mesoNetworkForceUpdate, Ftol, dPtol, P0, dt0);
	meso2Dobj.printMesoShearConfigCTCS2D(0.0);

	// get fixed contact network for G and B computation
	int NVTOT = meso2Dobj.getNVTOT();
	int NVVCTS = 0.5*NVTOT*(NVTOT-1);
	vector<bool> gijtmp(NVVCTS,0);

	// initialize temporary contact network
	meso2Dobj.getMesoVVContactNetwork(gijtmp);


	// double sxx1=meso2Dobj.getstress(0), syy1=meso2Dobj.getstress(1), sxy1=meso2Dobj.getstress(2);
	// meso2Dobj.mesoShearStrainEnthalpyFIRE(0.0, Ftol, P0, dt0, gijtmp);
	// // meso2Dobj.mesoEnthalpyFIRE(&meso2D::mesoNetworkForceUpdate, Ftol, dPtol, P0, dt0);
	// double sxx2=meso2Dobj.getstress(0), syy2=meso2Dobj.getstress(1), sxy2=meso2Dobj.getstress(2);

	// cout << "** Init: Sxx = " << sxx1 << ", Syy = " << syy1 << ", Sxy = " << sxy1 << endl;
	// cout << "** Scnd: Sxx = " << sxx2 << ", Syy = " << syy2 << ", Sxy = " << sxy2 << endl;
	// return 0;

	// object for shear strain
	meso2D mesoSaveObj(inputf,seed);
	mesoSaveObj = meso2Dobj;

	// -- Compute G numerically AT FIXED PRESSURE
	double dgamma = 1e-9;
	double gamma = 0.0;
	int NGAMMA = 10;
	int k = 0;

	// save shear stress and potential energy
	vector<double> UList(NGAMMA+1,0.0);
	vector<double> pList(NGAMMA+1,0.0);
	vector<double> sxyList(NGAMMA+1,0.0);
	vector<double> VList(NGAMMA+1,0.0);

	// save initial shear stress and potential energy
	UList.at(0) 		= meso2Dobj.getU();
	pList.at(0) 		= meso2Dobj.getPinst();
	sxyList.at(0) 		= meso2Dobj.getSinst();

	// loop over shear strains gamma, relax using FIRE, compute change in Sxy
	for (k=0; k<NGAMMA; k++){
		// update gamma for this iteration
		gamma += dgamma;

		// relax at fixed SHEAR STRAIN, PRESSURE, CONTACT NETWORK
		meso2Dobj.mesoShearStrainEnthalpyFIRE(gamma, Ftol, P0, dt0, gijtmp);

		// save mechanical information
		UList.at(k+1) 		= meso2Dobj.getU();
		pList.at(k+1) 		= meso2Dobj.getPinst();
		sxyList.at(k+1) 	= meso2Dobj.getSinst();

		// print
		meso2Dobj.printMesoShearConfigCTCS2D(gamma);
	}



	// -- Compute B numerically AT FIXED 

	// variables
	double Ptmp = P0;
	double dgammaB = 1e-2;
	gamma = 0.0;
	k = 0;

	// load old configuration
	meso2Dobj = mesoSaveObj;
	VList.at(0) = meso2Dobj.getL(0) * meso2Dobj.getL(1);

	// loop over strains in PRESSURE, relax using fire, compute change in volume
	for (k=0; k<NGAMMA; k++){
		// update gamma for this iteration
		gamma += dgammaB;

		// update temporary pressure based on gamma
		Ptmp = P0*(1 + gamma);

		// relax at fixed PRESSURE + CONTACT NETWORK
		meso2Dobj.mesoShearStrainEnthalpyFIRE(0.0, Ftol, Ptmp, dt0, gijtmp);

		// save volume
		VList.at(k+1) = meso2Dobj.getL(0) * meso2Dobj.getL(1);
	}


	// compute shear moduli using first derivatives of shear strain
	vector<double> G_sxy_list(NGAMMA-1,0.0);
	for (k=0; k<NGAMMA-1; k++)
		G_sxy_list.at(k) = -0.5*(sxyList.at(k+2) - sxyList.at(k))/dgamma;

	// compute bulk moduli using both normal finite diff and inverted
	vector<double> B_list(NGAMMA-1,0.0);
	vector<double> B_inv_list(NGAMMA-1,0.0);
	gamma = 0.0;
	double dP, Vtmp, V1, V2, dV, Vavg;
	for (k=0; k<NGAMMA-1; k++){
		// pressure difference
		dP = P0*dgammaB;

		// get volume information
		V1 = VList.at(k);
		V2 = VList.at(k+1);
		dV = V2 - V1;
		Vavg = 0.5*(V1 + V2);

		// approximate B
		B_list.at(k) = -Vavg*(dP/dV);
	}

	// print
	cout << "G modulus:" << endl;
	gamma = 0.0;
	cout << setprecision(12);
	for (k=0; k<NGAMMA-1; k++){
		cout << "k = " << k << ", gamma = " << gamma << ";   U = " << UList.at(k) << ";   P = " << pList.at(k) << ";    sxy = " << sxyList.at(k) << ";   G_sxy = " << G_sxy_list.at(k) << endl;
		gamma += dgamma;
	}

	cout << "B modulus:" << endl;
	gamma = 0.0;
	for (k=0; k<NGAMMA-1; k++){
		cout << "k = " << k << ", gamma = " << gamma << ";   V = " << VList.at(k) << ";   P = " << P0*(1 + gamma) << ";   B = " << B_list.at(k) << endl;
		gamma += dgammaB;
	}


	// say goodbye
	cout << "\n** Finished mesoModuliTest.cpp, ending. " << endl;

	return 0;
}