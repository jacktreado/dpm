#ifndef TUMOR_H
#define TUMOR_H

/*

	HEADER FILE FOR TUMOR CLASS

		-- Inherits from DPM class 
		-- For collections of 2D motile tumor cells
		-- Incorporates division (using preallocation) and adhesion
		-- Has own printing functions
			** Doesn't differentiate cell vs vertex info, just has vertex info with cell id in row

	Jack Treado, 06/04/21

*/

#include "dpm.h"


class tumor : public dpm{
protected:

	// temporary number of cells
	int ntmp;

	// adhesion parameters
	double l1, l2;

	// motility parameters
	double v0, Dr, Ds;

public:

	// Constructors and Destructors
	tumor(int n, int seed) : dpm(n,seed) { ntmp=0; l1=0.0; l2=0.0; v0=0.0; Dr=0.0; Ds=0.0; pbc[0]=0; pbc[1]=0; };

	// setters
	void setl1(double val) { l1 = val; };
	void setl2(double val) { l2 = val; };
	void setv0(double val) { v0 = val; };
	void setDr(double val) { Dr = val; };
	void setDs(double val) { Ds = val; };


	// initialize single cell
	void initializeSingleTumorCell();


	// editing & updating
	void divide(int ci);


	// force updates
	void tumorInteractionForces2D();


	// integrators
	void tumorFIRE2D(double Ftol, double dt0);


	// protocols
	void growMonolayer(double Ftol, double dt0);
	void activeBrownianCrawlers(int NT, double dt0);

};

#endif