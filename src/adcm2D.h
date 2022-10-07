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
		* Add bonding mechanism
		* Write force update function for circulolines in 2D
			* In normal shape force, add pressure and shear stress calculation
			* Write a shape force for bonded system
				* This will also need a pressure and shear stress
			* Add variables for bond length deletion cutoff
		* Incorporate constant-density vertex remapping (birth and death of i based on li and li-1)
		* Debug, make sure energy conserved with stick circuloline potential
		* Add compression-as-initialization protocol with FIRE for purely repulsive circlolines
		* Add crawling model
		* Add dynamic osmotic pressure
		* Simulate EMT via volume / adhesion changes
*/

#include "dpm.h"

class adcm2D : public dpm{
protected:

	// array of surface tensions
	std::vector<double> st;

	// arrays for bonds
	std::vector<int> bonds;
	std::vector<int> new_bonds;
	std::vector<int> labels;		
	std::vector<double> tproj;		// tproj[gi] is projection of vert gi onto edge bonds[gi]

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

	// initialization
	void nphCompression(double P0, double Ptol, double Ftol);

	// setup potential function pointer
	void useRepulsiveForce() { pwFrc = &adcm2D::SRRepulsivePWForce; }; 
	void useAttractiveForce() { pwFrc = &adcm2D::SRAttractivePWForce; }; 
	void useBondedForce() { pwFrc = &adcm2D::checkBondPWForce; };

	void useIndependentShapeForce() { shpFrc = &adcm2D::adcm2DShapeForces; };
	void useBondedShapeForce() { shpFrc = &adcm2D::adcm2DBondedShapeForces; };

	// circuloline geometry
	double getVertexEdgeProjection(const int gv, const int ge);
	double edge2VertexDistance(const int gv, const int ge, double &hx, double &hy, double &tev);

	// force update
	// NOTE: shape force is called before interaction force because, if bonds are present,
	// you should compute bonded shape forces before you change the network
	void adcm2DForceUpdate() { U = 0.0; Pinst = 0.0; Sinst = 0.0; fill(F.begin(), F.end(), 0.0); (*this.*shpFrc)(); circuloLinePWForceUpdate();  stressUpdate(); };

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

	// bond dynamics
	void updateBonds();
	void deleteVertex(const int gk);

	// -- Simulation Functions

	// NVE dynamics check
	void nve(int NT, double dt0, double T0, bool printDynamics);

	// FIRE relaxation at fixed pressure
	void nphFIRE(double Ftol, double dPtol, double P0, double dt0, bool printCompression);

	// printing
	void printADCM2DConfiguration();
};

#endif