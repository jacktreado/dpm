#ifndef MESO2D_H
#define MESO2D_H

/*

	HEADER FILE FOR MESO CLASS

		-- Inherits from DPM class 
		-- For collections of mesophyll cells
		-- Incorporates shape aging, stretching, box size changes
		-- ONLY FOR 2D

	Jack Treado, 06/09/21

*/


#include "dpm.h"

// global constants
const double kbmax = 1e-2;
const double calAmax = 2.0;

class meso2D;
typedef void (meso2D::*meso2DMemFn)(void);

class meso2D : public dpm{
protected:
	// int for maximum number of vertices
	int NVMAX;

	// bending energy per vertex
	// NOTE: will need to add different Hessian computation
	std::vector<double> kbi;

	// vertex-vertex contact network
	std::vector<bool> gij;

	// vv contacts per cell (both cell-level and vertex level)
	std::vector<int> zc;
	std::vector<int> zv;

	// pressure
	// NOTE: force routines will compute the partial dUdL first, need to add \sum_i F_i * r_i 
	// if you need instantaneous pressure contribution to minimized pressure
	double Pinst; 

	// adhesion parameters
	double betaEff;			// probability to break contact
	double ctcdel;			// influence of contact dependent adhesion
	double ctch;			// breaking strength

	// aging during development
	double cL; 		// aging of excess perimeter
	double aL; 		// distribution of aging to either contacting (0) or void (1) sections of perimeter
	double cB; 		// aging of preferred angles
	double cKb; 	// aging of bending stiffness

	// meso specific output objects
	std::ofstream hessout;
	std::ofstream ctcout;
public:

	// constructor and destructor
	meso2D(std::string &inputFile, int seed);
	meso2D(int n, int seed) : dpm(n,seed) { betaEff=0.0; ctcdel=1.0; ctch=0.5; cL=0.0; aL=1.0; cB=0.0; cKb=0.0; zc.resize(n); NVMAX = n; };

	// overloaded operators
	void operator=(const meso2D &rhs);

	// File openers
	void openHessObject(std::string& str) {
		hessout.open(str.c_str());
		if (!hessout.is_open()) {
			std::cout << "	ERROR: hessout could not open " << str << "..." << std::endl;
			exit(1);
		}
		else
			std::cout << "** Opening hess file " << str << " ..." << std::endl;
	}

	void openCTCObject(std::string& str) {
		ctcout.open(str.c_str());
		if (!ctcout.is_open()) {
			std::cout << "	ERROR: ctcout could not open " << str << "..." << std::endl;
			exit(1);
		}
		else
			std::cout << "** Opening ctc file " << str << " ..." << std::endl;
	}

	// getters
	double getPinst() { return Pinst; };


	// setters
	void setNVMAX(int val) { NVMAX = val; };
	void setkbi(double val) { fill(kbi.begin(), kbi.end(), val); };
	void setbetaEff(double val) { betaEff = val; };
	void setctcdel(double val) { ctcdel = val; };
	void setctch(double val) { ctch = val; };
	void setcL(double val) { cL = val; };
	void setaL(double val) { aL = val; };
	void setcB(double val) { cB = val; };
	void setcKb(double val) { cKb = val; };

	// editors & updates
	double meanl0();
	double meancalA0();
	double meant0();
	double meankb();

	// initialize mesophyll system
	void initializeVertexContactNetwork();
	void initializeMesophyllCells(double dispersion, double calA0, double phi0, double Ftol, int n1);

	// mesophyll cell interactions
	void initializeMesophyllBondNetwork();
	void mesoRepulsiveVertexForces();
	void mesoShapeForces();
	void mesoShapeForces(double gamma);
	void mesoNetworkForceUpdate();
	void mesoNetworkForceUpdate(double gamma, std::vector<bool> &gijtmp);
	void mesoPinForceUpdate(std::vector<double>& xpin, double kcspring);

	// integrators
	void mesoFIRE(meso2DMemFn forceCall, double Ftol, double dt0);
	void mesoEnthalpyFIRE(meso2DMemFn forceCall, double Ftol, double dPtol, double P0, double dt0);
	void mesoShearStrainEnthalpyFIRE(double gamma, double Ftol, double P0, double dt0, std::vector<bool> &gijtmp);
	void mesoPinFIRE(std::vector<double> &xpin, double Ftol, double dt0, double kcspring);
	void mesoNetworkNVE(std::ofstream &enout, meso2DMemFn forceCall, double T, double dt0, int NT, int NPRINTSKIP);
	void mesoShearStrainNVE(std::ofstream &enout, double gamma, double T, double dt0, int NT, int NPRINTSKIP, std::vector<bool> &gijtmp);

	// protocols
	void mesoNetworkExtension(meso2DMemFn forceCall, double Ftol, double dt0, double delShrink, double dphiPrint, double phiMin);
	void mesoPinExtension(double Ftol, double dt0, double hmax, double dh, double dhprint, double kcspring, int cellskip);
	void mesoFreeGrowth(meso2DMemFn forceCall, double Ftol, double dt0, double dl0, double da0, double dphiPrint, double a0max);
	void mesoNetworkEnthalpyMin(meso2DMemFn forceCall, double Ftol, double dPtol, double dt0, double da0, double dl0, double P0, double phiMin, int NMINSKIP);

	// protocol helpers
	void updateMesophyllBondNetwork(int CTCMIN, int PAIRMIN);
	void ageMesophyllShapeParameters();
	void relaxByAdding();
	void addMesophyllCellMaterial(double dl0);
	int mesoBondedCTCS(int gi);
	int mesoBondedPAIRS(int ci, int cj);
	void addVertex(int gi, double newl0);
	void t0ToCurrent();
	void t0ToReg();
	void getMesoVVContactNetwork(std::vector<bool> &gijtmp);
	double mesoInstantaneousPressure(std::vector<bool> &gijtmp);

	// hessian computation & linear response
	void mesoBendingHessian(Eigen::MatrixXd &Hb, Eigen::MatrixXd &Sb);
	void mesoSpringNetworkHessian(Eigen::MatrixXd &Hs, Eigen::MatrixXd &Ss);
	void mesoDynamicalMatrix(Eigen::MatrixXd &M, Eigen::MatrixXd &H, Eigen::MatrixXd &S);
	void mesoPrintLinearResponse(meso2DMemFn forceCall, double Ftol, double P0, double dt0);
	double numericalShearModulus(meso2DMemFn forceCall, double Ftol, double P0, double dt0);
	double numericalBulkModulus(meso2DMemFn forceCall, double Ftol, double P0, double dt0);

	// printing functions
	void printMesoNetwork2D();
	void printMesoNetworkCTCS2D();
	void printMesoPin2D(std::vector<double> &xpin, double h);
	void printMesoBondNetwork();
	void printMesoShearConfigCTCS2D(double gamma);
};

#endif