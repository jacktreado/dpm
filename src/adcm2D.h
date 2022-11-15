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
	double gam0;
	std::vector<double> st;

	// array of pressure differentials at each segment
	double Posm;
	std::vector<double> dp; 

	// adhesion ("adhering" circulolines decrease their surface tension)
	double W;

	// matrix of surface tensions between each cell
	std::vector< std::vector<double> > stMat;

	// // linking bonds arrays
	// std::vector<double *> link_loc;	 	// store locations of each linker on cell surface 
	// std::vector<int *> link_end; 		// store vertices that the ends of each linker are attached to, will need to update in add / delete 

	// pressure and shear stress
	double Pinst;
	double Sinst; 

	// pairwise interaction force function
	void (adcm2D::*pwFrc)(const int gi, const int gj);
	void (adcm2D::*shpFrc)();
public:

	// constructor and destructor
	adcm2D(int numcells, int numverts, double sizeDisp, double phi0, double boxLengthScale, double clScale, double initGamma0, int seed);
	adcm2D(std::string &inputFileStr, int seed, double boxLengthScale);

	// getters
	double getPinst() { return Pinst; };
	double getSinst() { return Sinst; };

	// setters
	void setgam0(double val) { gam0 = val; };
	void setgam(double val) { fill(st.begin(), st.end(), val); };
	void setstMat(double val) { for(int ci=0; ci<NCELLS+1; ci++) { fill(stMat.at(ci).begin(), stMat.at(ci).end(), val); } };
	void setW(double val) { W = val; };
	void setPosm(double val) { Posm = val; };

	// setup potential function pointer
	void useRepulsiveForce() { pwFrc = &adcm2D::SRRepulsivePWForce; }; 
	void useAttractiveForce() { pwFrc = &adcm2D::SRAttractivePWForce; }; 
	void useAlignmentForce() { pwFrc = &adcm2D::segmentPWAlignment; };
	void useBondedForce() { pwFrc = &adcm2D::checkBondPWForce; };

	void useIndependentShapeForce() { shpFrc = &adcm2D::adcm2DShapeForces; };
	void useContractileShapeForce() { shpFrc = &adcm2D::adcm2DContractileForces; };
	void useOsmoticPressureShapeForce() { shpFrc = &adcm2D::adcm2DOsmoticPressureShapeForces; };
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
	void SRRepulsivePWForce(const int gi, const int gj, bool&, bool&, bool&);
	void SRRepulsivePWForce(const int gi, const int gj) { bool a, b, c; SRRepulsivePWForce(gi, gj, a, b, c); };
	void SRAttractivePWForce(const int gi, const int gj, bool&, bool&, bool&);
	void SRAttractivePWForce(const int gi, const int gj) { bool a, b, c; SRAttractivePWForce(gi, gj, a, b, c); };
	void GhostAttractivePWForce(const int gi, const int gj);
	void segmentPWAlignment(const int gi, const int gj); 
	void alignmentForce(const int gi, const int gj);
	void checkBondPWForce(const int gi, const int gj);

	// shape forces
	void adcm2DShapeForces();
	void adcm2DContractileForces();
	void adcm2DOsmoticPressureShapeForces();
	void adcm2DBondedShapeForces();

	// stresses
	void stressUpdate();

	// -- Simulation Functions

	// NVE dynamics check
	void nve(int NT, double dt0, double T0, bool printDynamics);

	// FIRE relaxation at fixed pressure
	void nphFIRE(double Ftol, double dPtol, double P0, double dt0, bool printCompression);
	void nphSplitFIRE(double Ftol, double dPtol, double P0, double dt0, bool printCompression);
	void activeBrownianCrawling(const double Tsim, const double Tprint, const double dt0, const double v0, const double Dr, const double Ds);
	void activeTensionFluctuations(const double Tsim, const double Tprint, const double dt0, const double kneigh, const double deltaST);

	// printing
	void printADCM2DConfiguration();
};

#endif