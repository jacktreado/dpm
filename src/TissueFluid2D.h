#ifndef TISSUEFLUID2D_H
#define TISSUEFLUID2D_H

/*

	HEADER FILE FOR TissueFluid2D CLASS

	* Based on active deformable cells
    * Incorporates fluid flow between cells
    * Tracks fluid element label (-1 = res, >= 0 labels elements)
    
    TO-DO
        [ ] Add simple functions to traverse fluid_element_labels to compute things like element area, 
            contributions to A & B matrices for non-linear solve

        [ ] Write function to add fluid elements based on collision between two cells

        [ ] Add checks for merging / deleting elements based on cell separation

        [ ] Write function for non-linear solver, to determine Lagrange Multiplers (scaled pressures)
            to enforce area constraints

        [ ] Write test simulation for 2 cells to form channel with adhesion between them

    Issues / Anxieties
    
        * How to decide when to form / merge element? What if channel disappears, how to merge elements on either end?

	Jack Treado, 02/09/23

*/

#include "adcm2D.h"

class TissueFluid2D : public adcm2D {
protected:
    // fluid element label
    std::vector<int> fluid_element_label;

    // track pressure in each element
    std::vector<double> fluid_pressure;

public:

    // constructor: read in from file
    TissueFluid2D(std::string &inputFileStr, int seed, double boxLengthScale);

    // -- Simulation functions with constraints on single cell areas
	double areaConstraintLM(const int ci);
	void updateConstrainedAreas();
	void surfaceTensionForceUpdate();
};

#endif


