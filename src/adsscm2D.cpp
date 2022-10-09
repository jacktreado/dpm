#include "adsscm2D.h"

// namespacess
using namespace std;


adsscm2D::adsscm2D(int numcells, int numverts, double sizeDisp, double phi0, double boxLengthScale, double clScale, double gam0, int seed) : adcm2D(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, gam0, seed) {
	// initialize variables for adcm2D
	useRepulsiveForce();
	useIndependentShapeForce();
	v0 = 0.0; Dr = 1.0; Ds = 1.0;

	// initialize particles with gaussian sizes
	gaussian2D(sizeDisp, 1.0, numverts);
	st.resize(NVTOT);
	fill(st.begin(), st.end(), gam0);

	// initialize bond arrays
	bonds.resize(NVTOT);
	labels.resize(NVTOT);
	new_bonds.resize(NVTOT);
	tproj.resize(NVTOT);

	fill(bonds.begin(), bonds.end(), -1);
	fill(labels.begin(), labels.end(), -1);
	fill(new_bonds.begin(), new_bonds.end(), -1);
	fill(tproj.begin(), tproj.end(), 0.0);

	// rescale radii to be some fraction of l0
	for (int gi=0; gi<NVTOT; gi++)
		r.at(gi) *= clScale;

	// initialize spring constants
	setka(1.0);
	setkl(0.0);
	setkb(0.0);
	setkc(0.1);
	setl1(0.0);
	setl2(0.0);

	// initialize cell centers
	initializePositions2D(phi0, 1e-10);

	// initialize neighbor linked list
	initializeNeighborLinkedList2D(boxLengthScale);
}





// -- Shared Segment Forces

// pairwise force between two circulolines that may or may not be bonded
void adsscm2D::checkBondPWForce(const int gi, const int gj){
	// local variables
	bool vert_gi_on_edge_gj, vert_gj_on_edge_gi = 0;
	double hx_1, hy_1, hx_2, hy_2, hr_1, hr_2, t_1, t_2;

	// if not bonded, use sr repulsive pw force. 
	if (bonds.at(gi) != bonds.at(gj)){
		// update force
		SRRepulsivePWForce(gi, gj, vert_gi_on_edge_gj, vert_gj_on_edge_gi);

		// check whether to add to new_bonds list for later (vertex gi onto edge gj)
		if (vert_gi_on_edge_gj && bonds.at(gj) == -1 && bonds.at(gj) == -1){
			if (new_bonds.at(gi) == -1)
				new_bonds.at(gi) = gj;
			else {
				hr_1 = edge2VertexDistance(gi, gj, hx_1, hy_1, t_1);
				hr_2 = edge2VertexDistance(gi, new_bonds.at(gi), hx_2, hy_2, t_2);
				if (hr_1 < hr_2)
					new_bonds.at(gi) = gj;
			}
			
		}

		// same as above, but also check vertex gj onto edge gi
		if (vert_gj_on_edge_gi && bonds.at(gj) == -1 && bonds.at(gj) == -1){
			if (new_bonds.at(gj) == -1)
				new_bonds.at(gj) = gi;
			else {
				hr_1 = edge2VertexDistance(gj, gi, hx_1, hy_1, t_1);
				hr_2 = edge2VertexDistance(gj, new_bonds.at(gj), hx_2, hy_2, t_2);
				if (hr_1 < hr_2)
					new_bonds.at(gj) = gi;
			}
		}
	}
}

// shape forces for systems with shared segment dynamics
// TO DO
// -- copy from adcm2D::adcm2DShapeForces(), but include tij computation, shared forces



// update bonds based on new_bonds
void adsscm2D::updateBonds(){
	// local variables
	int gi;

	// loop over bonds, add new bonds, delete bonds from small segments
	for (gi=0; gi<NVTOT; gi++){
		// check to add bond at gi
		if (new_bonds.at(gi) > -1 && bonds.at(gi) == -1)
			bonds.at(gi) = new_bonds.at(gi);
	}

	// reset new_bonds to -1
	fill(new_bonds.begin(), new_bonds.end(), -1);
}

