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


class meso2D : public dpm{
protected:
	// bending energy per vertex
	// NOTE: will need to add different Hessian computation
	std::vector<double> kbi;

	// vertex-vertex contact network
	std::vector<bool> gij;

	// vv contacts per cell
	std::vector<int> z;

	// adhesion strength
	double betaEff;

	// aging during development
	double cL; 		// aging of excess perimeter
	double cB; 		// aging of preferred angles
	double cKb; 	// aging of bending stiffness
public:

	// constructor and destructor
	meso2D(int n, int seed) : dpm(n,seed) { betaEff=0.0; cL=0.0; cB=0.0; cKb=0.0; z.resize(n); };

	// setters
	void setkbi(double val) { fill(kbi.begin(), kbi.end(), val); };
	void setbetaEff(double val) { betaEff = val; };
	void setcL(double val) { cL = val; };
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
	void mesoShapeForces2D();
	void mesoNetworkForceUpdate();
	void mesoPinForceUpdate(std::vector<double>& xpin, double kcspring);

	// integrators
	void mesoNetworkFIRE(double Ftol, double dt0);
	void mesoPinFIRE(std::vector<double>& xpin, double Ftol, double dt0, double kcspring);
	void mesoNetworkNVE(double T, double dt0, int NT, int NPRINTSKIP);

	// protocols
	void mesoNetworkExtension(double Ftol, double dt0, double delShrink, double dphiPrint, double phiMin);
	void updateMesophyllBondNetwork();
	void ageMesophyllShapeParameters();
	void mesophyllPinExtension(double Ftol, double dt0, double hmax, double dh, double dhprint, double kcspring);


	// printing functions
	void printMesoNetwork2D();
	void printMesoPin2D(std::vector<double>& xpin, double h);
};

#endif