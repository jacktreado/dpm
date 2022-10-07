// test main file for meso shear modulus

// Compilation command:
// g++ -O3 --std=c++11 -I src main/test/mesoPressureTest.cpp src/*.cpp -o test.o

// header files
#include "meso2D.h"

// preprocessor macros
#define NDIM 2

using namespace std;

int main(){
	// local variables
	int seed = 1, NMINSKIP = 1;
	double Ftol = 1e-12, Ptol = 1e-6, dt0 = 0.005, dphi0 = 0.01;
	double boxLengthScale = 2.5, betaEff = 200.0, ctcdel = 1.0, ctch = 1, cL = 0.5, aL = 1.0, cB = 4.0, cKb = 0.0, kl = 1.0, kc = 1.0, kb0 = 0.4;
	double P0 = 1e-6, dPtol = 1e-10, phiMin = 0.9, da0 = 0.5, dl0 = 5.0, t0_min = 0.3;
	t0_min *= -0.5*PI;

	// pointer to dpm member function (should pt to null)
	dpmMemFn jammingForceUpdate = nullptr;
	meso2DMemFn mesoForceUpdate = nullptr;

	string inputf = "meso_n16.input";
	string posf = "pos.test";

	// instantiate object
	meso2D meso2Dobj(inputf,seed);

	// initialize neighbor linked list
	meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

	// open config, hessian and shear file
	meso2Dobj.openPosObject(posf);

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
	meso2Dobj.initializeMesophyllBondNetwork();
	meso2Dobj.t0ToCurrent();
	meso2Dobj.mesoFIRE(&meso2D::mesoNetworkForceUpdate, Ftol, dt0);
	meso2Dobj.t0ToCurrent();

	// run stretching simulation to create network
	// set max # of vertices
	int NVMAX = 5*meso2Dobj.getNVTOT();
	meso2Dobj.setNVMAX(NVMAX);
	meso2Dobj.mesoNetworkEnthalpyMin(&meso2D::mesoNetworkForceUpdate, Ftol, dPtol, dt0, da0, dl0, t0_min, P0, phiMin, NMINSKIP);

	// get fixed contact network for G and B computation
	int NVTOT = meso2Dobj.getNVTOT();
	int NVVCTS = 0.5*NVTOT*(NVTOT-1);
	vector<bool> gijtmp(NVVCTS,0);

	// initialize temporary contact network
	meso2Dobj.getMesoVVContactNetwork(gijtmp);

	// Loop over box sizes, measure pressure and total potential energy
	double L0 = meso2Dobj.getL(0);
	double dL = 1e-8*L0;
	int NCOMP = 20;
	double U, P, Pv, A;
	vector<double> UList(NCOMP+1,0.0);
	vector<double> PList(NCOMP+1,0.0);
	vector<double> PvList(NCOMP+1,0.0);
	vector<double> AList(NCOMP+1,0.0);

	// initialize lists
	AList.at(0) = meso2Dobj.getL(0) * meso2Dobj.getL(1);
	PvList.at(0) = 0.5*(meso2Dobj.getstress(0) + meso2Dobj.getstress(1))/AList.at(0);
	// P = meso2Dobj.mesoInstantaneousPressure(gijtmp);
	// U = meso2Dobj.getU();
	UList.at(0) = meso2Dobj.getU();
	PList.at(0) = meso2Dobj.getPinst();

	// loop
	for (int kk = 0; kk<NCOMP; kk++){
		// decrease box size
		meso2Dobj.setL(0,L0 - (kk+1)*dL);
		meso2Dobj.setL(1,L0 - (kk+1)*dL);
		meso2Dobj.initializeNeighborLinkedList2D(boxLengthScale);

		// relax potential energy at fixed volume
		meso2Dobj.mesoFIRE(&meso2D::mesoNetworkForceUpdate, Ftol, dt0);

		// save info
		A = meso2Dobj.getL(0) * meso2Dobj.getL(1);
		Pv = 0.5*(meso2Dobj.getstress(0) + meso2Dobj.getstress(1))/A;
		meso2Dobj.getMesoVVContactNetwork(gijtmp);
		// P = meso2Dobj.mesoInstantaneousPressure(gijtmp);
		// U = meso2Dobj.getU();
		U = meso2Dobj.getU();
		P = meso2Dobj.getPinst();

		// print to console for MATLAB
		cout << " ** kk = " << kk << ";  U = " << U << ";  P = " << P << ";  Pv = " << Pv << endl;

		// save
		UList.at(kk+1) = U;
		PList.at(kk+1) = P;
		PvList.at(kk+1) = Pv;
		AList.at(kk+1) = A;
	}

	// print to console
	for (int kk=0; kk<NCOMP+1; kk++)
		cout << kk << "  " << setprecision(15) << AList.at(kk) << "  " << UList.at(kk) << "  " << PList.at(kk) << "  " << PvList.at(kk) << endl;


	// say goodbye
	cout << "\n** Finished mesoPressureTest.cpp, ending. " << endl;

	return 0;
}







































