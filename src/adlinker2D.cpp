/*

	FUNCTION DEFINITIONS for ADLINKER2D class

    * Stands for Active Deformable Cell Model w/ linkers
	* only for use in two dimensions

	Jack Treado, 11/22/22

*/

#include "adlinker2D.h"

// namespaces
using namespace std;



/******************************

	C O N S T R U C T O R S,

	D E S T R U C T O R S,
	
	I N I T I A L I Z A T I O N

*******************************/

// constructor for adlinkers with normal initialization
adlinker2D::adlinker2D(int numcells, int numverts, double sizeDisp, double phi0, double boxLengthScale, double clScale, double initGamma0, int seed) : adcm2D(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, initGamma0, seed) {
	
}

// constructor for adlinkers with input initialization
adlinker2D::adlinker2D(std::string &inputFileStr, int seed, double boxLengthScale) : adcm2D(inputFileStr, seed, boxLengthScale) {
	
}



