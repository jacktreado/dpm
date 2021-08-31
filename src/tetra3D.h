#ifndef TETRA3D_H
#define TETRA3D_H


/*

	HEADER FILE FOR TETRA 3D CLASS

		-- Stores vertex information for NCELLS tetrahedra with n_\mu vertices each
		-- Primarily for 3D gelation simulations, could also pack to jamming
		-- Prints data to files


		NOTE that althought this class inherits from dpm class, many 2D-specific functions
		live in dpm class

		Cannibalizes the dynamical variables, but different structure for shape parameters
		and storage of vertices

	Jack Treado, 08/10/21

*/


#include "dpm.h"

class tetra3D : public dpm{
protected:

	// total number of faces and edges
	int NFTOT;
	int NETOT; 

	// face info
	std::vector<int> nf;		// # of faces per particle
	std::vector<int> finfo;		// mesh topology

	// preferred volume
	std::vector<int> v0;

public:

	

}