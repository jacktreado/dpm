#ifndef ADLINKER2D_H
#define ADLINKER2D_H

/*

	HEADER FILE FOR ADLINKER2D CLASS

	Active Deformable Cells with Linkers

		-- Inherits from ADCM2D class 
		-- For collections of active deformable cells * With Active Linker-based Adhesion * 
		-- only for use in 2D

	Jack Treado, 11/22/22 (HBD Mom!)
* 
	TO DO:
		* Simulate EMT via volume / adhesion changes
*/

#include "adcm2D.h"

class adlinker2D : public adcm2D {
protected:
	// linking bonds arrays
	std::vector<double *> link_loc;	 	// store locations of each linker on cell surface 
	std::vector<int *> link_end; 		// store vertices that the ends of each linker are attached to, will need to update in add / delete 
    std::vector<double> link_time;      // linker lifetime

    // time scales
    double t_on;                        // time to bind
    double t_off;                       // time to unbind
    double d_max;                       // maximum stretch distance

public:
    // constructor
    adlinker2D(int numcells, int numverts, double sizeDisp, double phi0, double boxLengthScale, double clScale, double initGamma0, int seed);
    adlinker2D(std::string &inputFileStr, int seed, double boxLengthScale);

    // add and subtract linkers
    void addLinkers();
    void delLinkers();

    // simulation
    void updateLinkerPositions();

    // print
    void printADLinker2DConfiguration();
};

#endif