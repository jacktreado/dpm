#ifndef ADCM2D_H
#define ADCM2D_H

/*

	HEADER FILE FOR ADCM2D CLASS

	Active Deformable Cell Model 

		-- Inherits from DPM class 
		-- For collections of active deformable cells * With Frictionless Boundaries *
		-- only for use in 2D

	Jack Treado, 09/13/22
* 
	TO DO:
		* Add bonding mechanism (NOW USING LINKING-BONDS. FOR SHARED SEGMENTS, SEE ADSSCM2D CLASS)
		* Write force update function for circulolines in 2D
			* In normal shape force, add shear stress calculation
		* Incorporate constant-density vertex remapping (birth and death of i based on li and li-1)
		* Debug, make sure energy conserved with stick circuloline potential
		* Add compression-as-initialization protocol with FIRE for purely repulsive circlolines
		* Add crawling model
		* Add dynamic osmotic pressure
		* Simulate EMT via volume / adhesion changes
*/

#include "dpm.h"

class adcm2D : public dpm {
protected:

	// preferred surface density of vertices
	double targetLength;

	// array of surface tensions
	std::vector<double> st;

	// // linking bonds arrays
	// std::vector<double *> link_loc;	 	// store locations of each linker on cell surface
	// std::vector<int *> link_end; 		// store vertices that the ends of each linker are attached to, will need to update in add / delete 

	// motility parameters
	double v0, Dr, Ds;
	std::vector<double> psi;

	// pressure and shear stress
	double Pinst;
	double Sinst; 

	// pairwise interaction force function
	void (adcm2D::*pwFrc)(const int gi, const int gj);
	void (adcm2D::*shpFrc)();

public:

	// constructor and destructor
	adcm2D(int numcells, int numverts, double sizeDisp, double phi0, double boxLengthScale, double clScale, double gam0, int seed);

	// setup potential function pointer
	void useRepulsiveForce() { pwFrc = &adcm2D::SRRepulsivePWForce; }; 
	void useAttractiveForce() { pwFrc = &adcm2D::SRAttractivePWForce; }; 
	void useBondedForce() { pwFrc = &adcm2D::checkBondPWForce; };

	void useIndependentShapeForce() { shpFrc = &adcm2D::adcm2DShapeForces; };
	void useBondedShapeForce() { shpFrc = &adcm2D::adcm2DBondedShapeForces; };

	// circuloline geometry
	double getVertexEdgeProjection(const int gv, const int ge);
	double edge2VertexDistance(const int gv, const int ge, double &hx, double &hy, double &tev);

	// dynamic cell vertices
	void checkVertices();
	void addVertex(const int gk);
	void deleteVertex(const int gk);

	// force update
	void adcm2DForceUpdate() { checkVertices(); U = 0.0; Pinst = 0.0; Sinst = 0.0; fill(F.begin(), F.end(), 0.0); (*this.*shpFrc)(); circuloLinePWForceUpdate(); stressUpdate(); };
	// void adcm2DForceUpdate() { U = 0.0; Pinst = 0.0; Sinst = 0.0; fill(F.begin(), F.end(), 0.0); (*this.*shpFrc)(); circuloLinePWForceUpdate(); stressUpdate(); };

	// interaction forces
	void circuloLinePWForceUpdate();
	void SRRepulsivePWForce(const int gi, const int gj, bool&, bool&);
	void SRRepulsivePWForce(const int gi, const int gj) { bool a, b; SRRepulsivePWForce(gi, gj, a, b); };
	void SRAttractivePWForce(const int gi, const int gj);
	void checkBondPWForce(const int gi, const int gj);

	// shape forces (ASSUMING kl = 0)
	void adcm2DShapeForces();
	void adcm2DBondedShapeForces();

	// stresses
	void stressUpdate();

	// -- Simulation Functions

	// NVE dynamics check
	void nve(int NT, double dt0, double T0, bool printDynamics);

	// FIRE relaxation at fixed pressure
	void nphFIRE(double Ftol, double dPtol, double P0, double dt0, bool printCompression);

	// printing
	void printADCM2DConfiguration();
};

#endif