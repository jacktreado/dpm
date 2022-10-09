#ifndef ADSSCM2D_H
#define ADSSCM2D_H

#include "adcm2D.h"

// Inherits from adcm2D, but just for shared segment dynamics

class adsscm2D : public adcm2D{
protected:

    // arrays for bonds
	std::vector<int> bonds;
	std::vector<int> new_bonds;
	std::vector<int> labels;		
	std::vector<double> tproj;		// tproj[gi] is projection of vert gi onto edge bonds[gi]

public:

adsscm2D(int numcells, int numverts, double sizeDisp, double phi0, double boxLengthScale, double clScale, double gam0, int seed);

void checkBondPWForce(const int gi, const int gj);
void updateBonds();

};

#endif