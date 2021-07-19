#ifndef TUMOR2D_H
#define TUMOR2D_H

/*

	HEADER FILE FOR TUMOR2D CLASS

		-- Inherits from DPM class 
		-- For collections of 2D motile tumor cells
		-- Incorporates division (using preallocation) and adhesion
		-- Has own printing functions
			** Doesn't differentiate cell vs vertex info, just has vertex info with cell id in row

	Jack Treado, 06/04/21

*/

#include "dpm.h"

class tumor2D;
typedef void (tumor2D::*tumor2DMemFn)(void);

class tumor2D : public dpm{
protected:

	// temporary number of tumor cells
	int tN;

	// surface tension info
	double gamtt;
	std::vector<double> l0_init;

	// wall pressures
	std::vector<double> wpress;

	// motility parameters
	double v0, Dr0, Ds;
	std::vector<double> psi;
	std::vector<double> Dr;

	// parameters to pin adipocytes
	double kecm, ecmbreak;
	std::vector<double> pinpos;
	std::vector<bool> pinattach;

public:

	// Constructors and Destructors
	tumor2D(std::string &inputFileString, int seed);
	tumor2D(int n, int seed) : dpm(n,seed) { tN=n; gamtt=0.0; v0=0.0; Dr0=0.0; Ds=0.0; kecm = 0.0; ecmbreak = 0.0; pbc[0]=0; pbc[1]=0; wpress.resize(2); }; 					// tumor-only constructor
	tumor2D(int n, int tNval, int seed) : dpm(n,seed) { tN=tNval; gamtt=0.0; v0=0.0; Dr0=0.0; Ds=0.0; kecm = 0.0; ecmbreak = 0.0; pbc[0]=0; pbc[1]=1; wpress.resize(2); }; 	// tumor+adipocyte constructor (x fixed, y periodic)


	// setters
	void setgamtt(double val) { gamtt = val; };
	void setv0(double val) { v0 = val; };
	void setDr0(double val) { Dr0 = val; };
	void setDs(double val) { Ds = val; };
	void setkecm(double val) { kecm = val; };
	void setecmbreak(double val) { ecmbreak = val; };
	void setl0_init();

	// initialization
	void initializeSingleTumorCell();
	void initializeTumorMonolayerPositions(double phi0, double Ftol, double kwell);
	void initializeTumorInterface(double aCalA0, double tCalA0, double aDisp, double tDisp, double areaRatio, int aNV, int tNV);
	void initializeTumorInterfacePositions(double phi0, double Ftol, double prt);

	// editing & updating
	void divide(int ci);

	// biology functions
	void psiDiffusion();
	void activeBrownianCrawlerUpdate();
	void updateECMAttachments(bool attach);
	void adipocyteECMAdhesionForces();

	// force updates
	void resetForcesAndEnergy();

	void repulsiveTumorForces();
	void stickyTumorForces();
	void repulsiveTumorInterfaceForces();
	void stickyTumorInterfaceForces();

	void repulsiveTumorForceUpdate();
	void stickyTumorForceUpdate();
	void repulsiveTumorInterfaceForceUpdate();
	void stickyTumorInterfaceForceUpdate();

	// integrators
	void tumorFIRE(tumor2DMemFn forceCall, double Ftol, double dt0);

	// protocols
	void setupCheck();
	void tumorCompression(double Ftol, double Ptol, double dt0, double dphi0);
	void invasion(tumor2DMemFn forceCall, double dDr, double dPsi, double Drmin, int NT, int NPRINTSKIP);
	void crawling(tumor2DMemFn forceCall, int NT, int NPRINTSKIP);

	// print functions
	void printTumorInterface(double t);
	void printTumorCells(double t);
};

#endif