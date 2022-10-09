/*

	FUNCTION DEFINITIONS for ADCM2D class

    * Stands for Active Deformable Cell Model 
	* only for use in two dimensions

	Jack Treado, 09/13/22

*/

#include "adcm2D.h"

// namespacess
using namespace std;



/******************************

	C O N S T R U C T O R S  & 
	
	I N I T I A L I Z A T I O N

*******************************/

// constructor for initializing at low density
// * numcells: # of cells in the simulation
// * numverts: # of vertices for average cell size
// * sig0: dispersion of radii / vertex number
// * phi0: initial packing fraction of seed disk particles particles
// * gam0: initial mean surface tension
// * seed: seed for randum number generator
adcm2D::adcm2D(int numcells, int numverts, double sizeDisp, double phi0, double boxLengthScale, double clScale, double gam0, int seed) : dpm(numcells, 2, seed) {
	// initialize variables for adcm2D
	useRepulsiveForce();
	useIndependentShapeForce();
	v0 = 0.0; Dr = 1.0; Ds = 1.0;

	// initialize particles with gaussian sizes
	gaussian2D(sizeDisp, 1.0, numverts);
	st.resize(NVTOT);
	fill(st.begin(), st.end(), gam0);

	// target length (based on circle with numverts equally-spaced & unit area)
	targetLength = (2.0 * sqrt(PI)) / numverts ;

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




/******************************

	C I R C U L O - L I N E

	G E O M E T R Y

*******************************/

double adcm2D::getVertexEdgeProjection(const int gv, const int ge){
	// local variables
	double dx, dy, ux, uy, l;

	// get vector pointing from edge to vertex
	dx = deltaX(ge, gv, 0);
	dy = deltaX(ge, gv, 1);

	// segment length
	l = seg(ge);

	// unit vector
	ux = segX(ge, 0) / l;
	uy = segX(ge, 1) / l;

	// projection
	return (dx * ux + dy * uy) / l;
}

double adcm2D::edge2VertexDistance(const int gv, const int ge, double &hx, double &hy, double &tev){
	// local variables
	double cp, dp, dx, dy, ux, uy, l;

	// get vector pointing from edge to vertex
	dx = deltaX(ge, gv, 0);
	dy = deltaX(ge, gv, 1);

	// unit vector along segment
	ux = unitSegX(ge, 0);
	uy = unitSegX(ge, 1);

	// cross product
	cp = dx * uy - dy * ux;

	// vector distance from edge to vertex
	hx = uy * cp;
	hy = -ux * cp;

	// projection
	tev = (dx * ux + dy * uy) / seg(ge);

	// return distance
	return sqrt(hx * hx + hy * hy);
}








/******************************

	C I R C U L O - L I N E

	F O R C E S

*******************************/


// force update for circulolines
void adcm2D::circuloLinePWForceUpdate(){
	// local variables
	int bi, bj, pi, pj, gi, gj, ci, vi, cj, vj;

	// sort particles
	sortNeighborLinkedList2D();

	// loop over cell collision linked list boxes
	for (bi = 0; bi < NBX; bi++) {
		// get start of list of vertices
		pi = head[bi];

		// loop over linked list
		while (pi > 0) {
			// real particle index
			gi = pi - 1;

			// cell index of gi
			cindices(ci, vi, gi);

			// next particle in list
			pj = list[pi];

			// loop down neighbors of pi in same cell
			while (pj > 0) {
				// real index of pj
				gj = pj - 1;

				// cell index of gj
				cindices(cj, vj, gj);

				if (ci == cj) {
					pj = list[pj];
					continue;
				}

				// compute force between circulolines
				// SRRepulsivePWForce(gi, gj);
				(*this.*pwFrc)(gi, gj);

				// update pj
				pj = list[pj];
			}

			// test overlaps with forward neighboring cells
			for (bj = 0; bj < NNN; bj++) {
				// only check if boundaries permit
				if (nn[bi][bj] == -1)
					continue;

				// get first particle in neighboring cell
				pj = head[nn[bi][bj]];

				// loop down neighbors of pi in same cell
				while (pj > 0) {
					// real index of pj
					gj = pj - 1;

					if (gj == ip1[gi] || gj == im1[gi]) {
						pj = list[pj];
						continue;
					}

					// compute force between circulolines
					// SRRepulsivePWForce(gi, gj);
					(*this.*pwFrc)(gi, gj);

					// update pj
					pj = list[pj];
				}
			}

			// update pi index to be next
			pi = list[pi];
		}
	}
}

// short-range (SR) repulsive pairwise force between two circulolines
// TO DO
// -- Debug segment to segment forces
void adcm2D::SRRepulsivePWForce(const int gi, const int gj, bool &vert_gi_on_edge_gj, bool &vert_gj_on_edge_gi){
	// local variables
	double tij, tji;
	double hr_i2j, hr_j2i, hx_i2j, hx_j2i, hy_i2j, hy_j2i;
	double ftmp, dfx, dfy;

	// contact distance
	double sij = r[gi] + r[gj];

	// indices
	int gip1 = ip1[gi];
	int gim1 = im1[gi];
	int gjp1 = ip1[gj];
	int gjm1 = im1[gj];

	// distance between vertex i and j
	double dx = deltaX(gi, gj, 0);
	double dy = deltaX(gi, gj, 1);
	double dr = sqrt(dx * dx + dy * dy);

	// Gut check distance
	double lij = 0.5 * (seg(gi) + seg(gj));
	double dlijApprox = 0.5 * (dr + deltaR(gip1, gjp1));
	if (dlijApprox < 80.0 * lij) {
		// projection of vertex j onto edge i
		hr_j2i = edge2VertexDistance(gj, gi, hx_j2i, hy_j2i, tij);

		// projection of vertex i onto edge j
		hr_i2j = edge2VertexDistance(gi, gj, hx_i2j, hy_i2j, tji);
		
		// check determine force based on projection
		if (tij < 0 && dr < sij){
			// force info
			ftmp = -(kc / sij) * (1 - (dr / sij));	// NOTE: -1 x mag of the force, it keeps to convention that dfx is > 0 for forces on gi
			dfx = ftmp * (dx / dr);
			dfy = ftmp * (dy / dr);

			// add to force
			addF(dfx, gi, 0);
			addF(dfy, gi, 1);
			addF(-dfx, gj, 0);
			addF(-dfy, gj, 1);

			// add to potential energy
			addU(0.5 * kc * pow(1 - (dr / sij), 2.0));

			// add to pressure
			Pinst += ftmp * (dr / L[0]);
		}
		else if (tij > 0 && tij < 1 && hr_j2i < sij) {
			// report
			vert_gj_on_edge_gi = 1;

			// force info
			ftmp = -(kc / sij) * (1 - (hr_j2i / sij));	// NOTE: this is -1 x mag of the force, it just sticks to convention that dfx is > 0 for forces on gi
			dfx = ftmp * (hx_j2i / hr_j2i);
			dfy = ftmp * (hy_j2i / hr_j2i);

			// add to force (distribute so torque is balanced)
			addF(dfx * (1 - tij), gi, 0);
			addF(dfy * (1 - tij), gi, 1);
			addF(dfx * tij, gip1, 0);
			addF(dfy * tij, gip1, 1);
			addF(-dfx, gj, 0);
			addF(-dfy, gj, 1);

			// add to potential energy
			addU(0.5 * kc * pow(1 - (hr_j2i / sij), 2.0));

			// add to pressure
			Pinst += ftmp * (hr_j2i / L[0]);
		}

		// also check projection from edge j to vertex i (same as tij case, but note use of -dx and -dy)
		if (tji > 0 && tji < 1 && hr_i2j < sij && false) {
			// report
			vert_gi_on_edge_gj = 1;

			// force info
			ftmp = -(kc / sij) * (1 - (hr_i2j / sij));	// NOTE: this is -1 x mag of the force, it just sticks to convention that dfx is > 0 for forces on gi
			dfx = ftmp * (hx_i2j / hr_i2j);
			dfy = ftmp * (hy_i2j / hr_i2j);

			// add to force (distribute so torque is balanced)
			addF(-dfx, gi, 0);
			addF(-dfy, gi, 1);
			addF(dfx * (1 - tji), gj, 0);
			addF(dfy * (1 - tji), gj, 1);
			addF(dfx * tji, gjp1, 0);
			addF(dfy * tji, gjp1, 1);
			
			// add to potential energy
			addU(0.5 * kc * pow(1 - (hr_i2j / sij), 2.0));

			// add to pressure
			Pinst += ftmp * (hr_i2j / L[0]);
		}
	}
}

/*
// short-range (SR) attractive pairwise force between two circulolines
void adcm2D::SRAttractivePWForce(const int gi, const int gj){
	// local variables
	double tij, tji;
	double hr_i2j, hr_j2i, hx_i2j, hx_j2i, hy_i2j, hy_j2i;
	double ftmp, dfx, dfy;

	// contact distance
	double sij = r[gi] + r[gj];

	// indices
	int gip1 = ip1[gi];
	int gim1 = im1[gi];
	int gjp1 = ip1[gj];
	int gjm1 = im1[gj];

	// distance between vertex i and j
	double dx = deltaX(gi, gj, 0);
	double dy = deltaX(gi, gj, 1);
	double dr = sqrt(dx * dx + dy * dy);

	// Gut check distance
	double lij = 0.5 * (seg(gi) + seg(gj));
	double dlijApprox = 0.5 * (dr + deltaR(gip1, gjp1));
	if (dlijApprox < 80.0 * lij) {
		// projection of vertex j onto edge i
		hr_j2i = edge2VertexDistance(gj, gi, hx_j2i, hy_j2i, tij);

		// projection of vertex i onto edge j
		hr_i2j = edge2VertexDistance(gi, gj, hx_i2j, hy_i2j, tji);
		
		// check determine force based on projection
		if (tij < 0 && dr < sij){
			// force info
			ftmp = -(kc / sij) * (1 - (dr / sij));	// NOTE: -1 x mag of the force, it keeps to convention that dfx is > 0 for forces on gi
			dfx = ftmp * (dx / dr);
			dfy = ftmp * (dy / dr);

			// add to force
			addF(dfx, gi, 0);
			addF(dfy, gi, 1);
			addF(-dfx, gj, 0);
			addF(-dfy, gj, 1);

			// add to potential energy
			addU(0.5 * kc * pow(1 - (dr / sij), 2.0));

			// add to pressure
			Pinst += ftmp * (dr / L[0]);
		}
		else if (tij > 0 && tij < 1 && hr_j2i < sij) {
			// force info
			ftmp = -(kc / sij) * (1 - (hr_j2i / sij));	// NOTE: this is -1 x mag of the force, it just sticks to convention that dfx is > 0 for forces on gi
			dfx = ftmp * (hx_j2i / hr_j2i);
			dfy = ftmp * (hy_j2i / hr_j2i);

			// add to force (distribute so torque is balanced)
			addF(dfx * (1 - tij), gi, 0);
			addF(dfy * (1 - tij), gi, 1);
			addF(dfx * tij, gip1, 0);
			addF(dfy * tij, gip1, 1);
			addF(-dfx, gj, 0);
			addF(-dfy, gj, 1);

			// add to potential energy
			addU(0.5 * kc * pow(1 - (hr_j2i / sij), 2.0));

			// add to pressure
			Pinst += ftmp * (hr_j2i / L[0]);
		}

		// also check projection from edge j to vertex i (same as tij case, but note use of -dx and -dy)
		if (tji > 0 && tji < 1 && hr_i2j < sij && false) {
			// force info
			ftmp = -(kc / sij) * (1 - (hr_i2j / sij));	// NOTE: this is -1 x mag of the force, it just sticks to convention that dfx is > 0 for forces on gi
			dfx = ftmp * (hx_i2j / hr_i2j);
			dfy = ftmp * (hy_i2j / hr_i2j);

			// add to force (distribute so torque is balanced)
			addF(-dfx, gi, 0);
			addF(-dfy, gi, 1);
			addF(dfx * (1 - tji), gj, 0);
			addF(dfy * (1 - tji), gj, 1);
			addF(dfx * tji, gjp1, 0);
			addF(dfy * tji, gjp1, 1);
			
			// add to potential energy
			addU(0.5 * kc * pow(1 - (hr_i2j / sij), 2.0));

			// add to pressure
			Pinst += ftmp * (hr_i2j / L[0]);
		}
	}
}
*/


// shape forces for systems with only surface tension (ASSUMING kb, kl == 0)
void adcm2D::adcm2DShapeForces(){
	// local variables
	int ci, gi, vi, nvtmp;
	double fa, cx, cy, xi, yi;
	double a0tmp, atmp;
	double dx, dy, da, dli, dlim1;
	double lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
	double rim1x, rim1y, rix, riy, rip1x, rip1y;

	// loop over vertices, add to force
	ci = 0;
	for (gi = 0; gi < NVTOT; gi++) {

		// -- Area force (and get cell index ci)
		if (ci < NCELLS) {
			if (gi == szList[ci]) {
				// shape information
				nvtmp = nv.at(ci);
				a0tmp = a0.at(ci);

				// compute area deviation
				atmp = area(ci);
				da = (atmp / a0tmp) - 1.0;

				// area force constant
				fa = (ka / a0tmp) * da;

				// update potential energy, pressure
				// NOTE: this is not yet the pressure, just the contribution from \partial U / \partial L
				U 		+= 0.5 * ka * (da * da);
				Pinst 	+= (2.0 * fa * atmp) / L[0];

				// compute cell center of mass
				xi = x[NDIM * gi];
				yi = x[NDIM * gi + 1];
				cx = xi;
				cy = yi;
				for (vi = 1; vi < nvtmp; vi++) {
					// get distances between vim1 and vi
					dx = x[NDIM * (gi + vi)] - xi;
					dy = x[NDIM * (gi + vi) + 1] - yi;
					if (pbc[0])
						dx -= L[0] * round(dx / L[0]);
					if (pbc[1])
						dy -= L[1] * round(dy / L[1]);

					// add to centers
					xi += dx;
					yi += dy;

					cx += xi;
					cy += yi;
				}
				cx /= nvtmp;
				cy /= nvtmp;

				// get coordinates relative to center of mass
				rix = x.at(NDIM * gi) - cx;
				riy = x.at(NDIM * gi + 1) - cy;

				rim1x = x.at(NDIM * im1[gi]) - cx;
				rim1y = x.at(NDIM * im1[gi] + 1) - cy;
				if (pbc[0])
					rim1x -= L[0] * round(rim1x / L[0]);
				if (pbc[1])
					rim1y -= L[1] * round(rim1y / L[1]);

				lim1x = rix - rim1x;
				lim1y = riy - rim1y;
				lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);

				// increment cell index
				ci++;
			}
		}

		// get next adjacent vertices
		rip1x = x[NDIM * ip1[gi]] - cx;
		rip1y = x[NDIM * ip1[gi] + 1] - cy;
		if (pbc[0])
			rip1x -= L[0] * round(rip1x / L[0]);
		if (pbc[1])
			rip1y -= L[1] * round(rip1y / L[1]);

		// -- Area force

		// add to force
		F[NDIM * gi] += 0.5 * fa * (rim1y - rip1y);
		F[NDIM * gi + 1] += 0.5 * fa * (rip1x - rim1x);

		// -- Surface tension force

		// segment lengths
		lix = rip1x - rix;
		liy = rip1y - riy;
		li = sqrt(lix * lix + liy * liy);

		// add to force
		F[NDIM * gi] += (st[gi] * (lix / li)) - (st[im1[gi]] * (lim1x / lim1));
		F[NDIM * gi + 1] += (st[gi] * (liy / li)) - (st[im1[gi]] * (lim1y / lim1));

		// add to potential energy
		U += st[gi] * li;

		// add to pressure and shear stress (TO DO: add shear stress)
		// NOTE: this is not yet the pressure, just the contribution from \partial U / \partial L
		Pinst += (st[gi] * li) / L[0];

		// update old coordinates
		rim1x = rix;
		rix = rip1x;

		rim1y = riy;
		riy = rip1y;

		// update segment vectors
		lim1x = lix;
		lim1y = liy;
		lim1 = li;
	}
}


// update stresses
void adcm2D::stressUpdate(){
	// local variables
	int i;
	double Fvdiag, Fvxy;

	// virial force contribution to pressure
	Fvdiag = 0.0;
	Fvxy = 0.0;
	for (i=0; i<vertDOF; i++){
		Fvdiag += F[i] * x[i];
		if (i % NDIM == 0)
			Fvxy += F[i] * x[i+1];
	}
	Fvdiag /= 2.0 * L[0] * L[1];
	Fvxy /= L[0] * L[1];

	// finalize computation of pressure
	Pinst = Fvdiag - (Pinst / (2.0 * L[0]));
	Sinst = Fvxy - (Sinst / (L[0] * L[1]));
}


/******************************

	B O N D  &  V E R T E X

	D Y N A M I C S

*******************************/

// check vertices to see if some should be deleted or added
void adcm2D::checkVertices(){
	// local variables
	int gi;
	double li;

	// lists for deletion and addition
	vector<int> verts2Del;
	vector<int> verts2Add;

	// loop over all vertices, check length, delete first
	for (gi=0; gi<NVTOT; gi++){
		// get segment length
		li = seg(gi);
		if (li < 0.5 * targetLength)
			verts2Del.push_back(gi);
	}

	// loop over vertices to delete, delete them
	for (gi=0; gi<verts2Del.size(); gi++)
		deleteVertex(verts2Del.at(gi) - gi);

	// loop over all vertices 
	for (gi=0; gi<NVTOT; gi++){
		// get new segment length
		li = seg(gi);
		if (li > 2.0 * targetLength)
			verts2Add.push_back(gi);
	}

	// loop over vertices to add, add them
	for (gi=0; gi<verts2Add.size(); gi++)
		addVertex(verts2Add.at(gi) + gi);
}

// delete vertex, reorganize indexing
void adcm2D::deleteVertex(const int gk){
	// local variables
	int ck, vk, ci, vi, vim1, vip1;
	int gi, d;

	// get cell level information
	cindices(ck, vk, gk);

	// loop over all vertices, reshuffle and overwrite, pop back all vectors
	for (gi=gk+1; gi<NVTOT; gi++){
		// shift things in NVTOT x 1 arrays
		r.at(gi - 1) = r.at(gi);
		st.at(gi - 1) = st.at(gi);

		// shift things in NVTOT x NDIM arrays
		for (d=0; d<NDIM; d++){
			x.at(NDIM * (gi - 1) + d) = x.at(NDIM * gi + d);
			v.at(NDIM * (gi - 1) + d) = v.at(NDIM * gi + d);
			F.at(NDIM * (gi - 1) + d) = F.at(NDIM * gi + d);
		}
	}

	// delete entries in vectors
	list.pop_back();
	r.pop_back();
	st.pop_back();
	im1.pop_back();
	ip1.pop_back();
	for (d=0; d<NDIM; d++){
		x.pop_back();
		v.pop_back();
		F.pop_back();
	}

	// change NVTOT
	NVTOT--;
	vertDOF -= NDIM;
	nv.at(ck)--;

	// recompute szList
	for (ci=1; ci<NCELLS; ci++)
		szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

	// resort ip1 and im1 based on new number of vertices
	for (ci=0; ci<NCELLS; ci++) {
		// vertex indexing
		for (vi = 0; vi < nv.at(ci); vi++) {
			// wrap local indices
			vim1 = (vi - 1 + nv.at(ci)) % nv.at(ci);
			vip1 = (vi + 1) % nv.at(ci);

			// get global wrapped indices
			gi = gindex(ci, vi);
			im1.at(gi) = gindex(ci, vim1);
			ip1.at(gi) = gindex(ci, vip1);
		}
	}
}

// add vertex, reorganize indexing
void adcm2D::addVertex(const int gk){
	// local variables
	int ck, vk, ci, vi, vim1, vip1;
	int gi, d;
	double dx, dy, xnew, ynew;

	// get cell level information
	cindices(ck, vk, gk);

	// determine new vertex position
	dx = segX(gk, 0);
	dy = segX(gk, 1);
	xnew = getx(gk, 0) + 0.5 * dx;
	ynew = getx(gk, 1) + 0.5 * dy;

	// add entries in vectors
	list.push_back(0);
	r.push_back(0.0);
	st.push_back(0.0);
	im1.push_back(0);
	ip1.push_back(0);
	for (d=0; d<NDIM; d++){
		x.push_back(0.0);
		v.push_back(0.0);
		F.push_back(0.0);
	}

	// loop over all vertices, reshuffle and overwrite, push back all vectors
	for (gi = NVTOT-1; gi > gk; gi--){
		// shift things in NVTOT x 1 arrays
		r.at(gi + 1) = r.at(gi);
		st.at(gi + 1) = st.at(gi);

		// shift things in NVTOT x NDIM arrays
		for (d=0; d<NDIM; d++){
			x.at(NDIM * (gi + 1) + d) = x.at(NDIM * gi + d);
			v.at(NDIM * (gi + 1) + d) = v.at(NDIM * gi + d);
			F.at(NDIM * (gi + 1) + d) = F.at(NDIM * gi + d);
		}
	}
		

	// change NVTOT
	NVTOT++;
	vertDOF += NDIM;
	nv.at(ck)++;

	// recompute szList
	for (ci = 1; ci < NCELLS; ci++)
		szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

	// resort ip1 and im1 based on new number of vertices
	for (ci = 0; ci < NCELLS; ci++) {
		// vertex indexing
		for (vi = 0; vi < nv.at(ci); vi++) {
			// wrap local indices
			vim1 = (vi - 1 + nv.at(ci)) % nv.at(ci);
			vip1 = (vi + 1) % nv.at(ci);

			// get global wrapped indices
			gi = gindex(ci, vi);
			im1.at(gi) = gindex(ci, vim1);
			ip1.at(gi) = gindex(ci, vip1);
		}
	}

	// new vertex position
	x[NDIM*(gk + 1)] = xnew;
	x[NDIM*(gk + 1) + 1] = ynew;

	// new forces + velocity
	v[NDIM*(gk + 1)] = 0.0;
	v[NDIM*(gk + 1) + 1] = 0.0;
	F[NDIM*(gk + 1)] = 0.0;
	F[NDIM*(gk + 1) + 1] = 0.0;

	// surface tension and radius
	r[gk + 1] = r[gk];
	st[gk + 1] = st[gk];
}










/******************************

	S I M U L A T I O N

	F U N C T I O N S

*******************************/


// NVE dynamics, report total energy during simulation 
void adcm2D::nve(int NT, double dt0, double T0, bool printDynamics){
	// local variables
	int tt, d;

	// set dt
	dt = dt0;

	// initialize energies
	adcm2DForceUpdate();
	drawVelocities2D(T0);

	// loop over time, integate VV
	for (tt=0; tt<NT; tt++){
		// VV velocity half-step update
		vvVelUpdate();

		// VV positions full-step update
		vvPosUpdate();

		// Update Forces & Potential Energy
		adcm2DForceUpdate();

		// VV velocity half-step update
		vvVelUpdate();

		// print statement
		if (tt % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << " 	N V E 									" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** tt 		= " << tt << endl;
			cout << "	** U 		= " << U << endl;
			cout << " 	** K 		= " << vertexKineticEnergy() << endl;
			cout << " 	** E 		= " << U + vertexKineticEnergy() << endl;
			cout << "	** dt 		= " << dt << endl;

			// print it
			if (printDynamics && posout.is_open())
				printADCM2DConfiguration();
		}
	}

}


// Enthalpy minimization
void adcm2D::nphFIRE(double Ftol, double dPtol, double P0, double dt0, bool printCompression){
	// local variables
	int i, d;

	// box size, momentum, and internal virial pressure
	double V=L[0]*L[1], Pi=0.0, P=0.0;

	// check to see if cell linked-list has been initialized
	if (NBX == -1){
		cerr << "	** ERROR: In adcm2D::mesoEnthalpyFIRE, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
		exit(1);
	}

	// FIRE variables
	double PFIRE, fnorm, fcheck, fcheckFrc, dPcheck, vnorm, alpha, dtmax, dtmin;
	int npPos, npNeg, fireit;

	// set dt based on geometric parameters
	dt = dt0;

	// Initialize FIRE variables
	PFIRE  		= 0;	
	fnorm 		= 0;
	vnorm 		= 0;
	alpha   	= alpha0;

	dtmax   	= 20.0 * dt;
	dtmin   	= 0.1 * dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;
	fcheckFrc 	= fcheck;

	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);
	adcm2DForceUpdate();
	P = Pinst;

	// relax forces using FIRE
	while ((fcheck > Ftol || dPcheck > dPtol || fireit < NDELAY) && fireit < itmax){
		// VV VELOCITY UPDATE #1
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*dt*(F[i] - 0.5*v[i]*(Pi/V));
		Pi += 0.5*dt*(P-P0);

		// compute PFIRE
		PFIRE = 0.0;
		for (i=0; i<vertDOF; i++)
			PFIRE += v[i]*(F[i] - 0.5*v[i]*(Pi/V));
		PFIRE += Pi*(P-P0);

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << " 	F I R E 								" << endl;
			cout << " 	E N T H A L P Y  						" << endl;
			cout << "	M I N I M I Z A T I O N 				" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit 	= " << fireit << endl;
			cout << "	** fcheck 	= " << fcheck << endl;
			cout << " 	** fcheckFrc 	= " << fcheckFrc << endl;
			cout << "	** U 		= " << U << endl;
			cout << "	** dt 		= " << dt << endl;
			cout << "	** PFIRE 	= " << PFIRE << endl;
			cout << "	** alpha 	= " << alpha << endl;
			cout << "	** npPos 	= " << npPos << endl;
			cout << "	** npNeg 	= " << npNeg << endl;
			cout << "	** V  		= " << V << endl;
			cout << "	** P 		= " << P << endl;
			cout << "	** Pi 		= " << Pi << endl;

			// print it
			if (printCompression && posout.is_open())
				printADCM2DConfiguration();
		}

		// Adjust simulation based on net motion of degrees of freedom
		if (PFIRE > 0){
			// increase positive counter
			npPos++;

			// reset negative counter
			npNeg = 0;

			// alter simulation if enough positive steps have been taken
			if (npPos > NDELAY){
				// change time step
				if (dt*finc < dtmax)
					dt *= finc;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else{
			// reset positive counter
			npPos = 0;

			// increase negative counter
			npNeg++;

			// check if simulation is stuck
			if (npNeg > NNEGMAX){
				cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
				exit(1);
			}

			// take half step backwards, reset velocities
			for (i=0; i<vertDOF; i++){
				// take half step backwards
				x[i] -= 0.5*dt*(v[i] + 0.5*x[i]*(Pi/V));

				// reset vertex velocities
				v[i] = 0.0;
			}
			V -= 0.5*dt*Pi;

			// decrease time step if past initial delay
			if (fireit > NDELAY){
				// decrease time step 
				if (dt*fdec > dtmin)
					dt *= fdec;

				// reset alpha
				alpha = alpha0;
			}
		}

		// compute fnorm, vnorm and P
		fnorm = 0.0;
		vnorm = 0.0;
		for (i=0; i<vertDOF; i++){
			fnorm 	+= pow(F[i] - 0.5*v[i]*(Pi/V),2.0);
			vnorm 	+= v[i]*v[i];
		}
		fnorm += pow(P - P0,2.0);
		vnorm += Pi*Pi;
		fnorm = sqrt(fnorm);
		vnorm = sqrt(vnorm);

		// update velocities (s.d. vs inertial dynamics) only if forces are acting
		if (fnorm > 0){
			for (i=0; i<vertDOF; i++)
				v[i] = (1 - alpha)*v[i] + alpha*((F[i] - 0.5*v[i]*(Pi/V))/fnorm)*vnorm;
			Pi = (1 - alpha)*Pi + alpha*((P - P0)/fnorm)*vnorm;
		}

		// VV POSITION UPDATE
		for (i=0; i<vertDOF; i++)
			x[i] += dt*(v[i] + 0.5*x[i]*(Pi/V));
		V += dt*Pi;
		L[0] = sqrt(V);
		L[1] = sqrt(V);
		for (d = 0; d < NDIM; d++)
			lb[d] = L[d] / sb[d];

		// update forces, pressure
		adcm2DForceUpdate();
		P = Pinst;

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*dt*(F[i] - 0.5*v[i]*(Pi/V));
		Pi += 0.5*dt*(P - P0);

		// update fcheck based on fnorm (= force per degree of freedom)
		fcheckFrc = 0.0;
		for (i=0; i<vertDOF; i++)
			fcheckFrc += pow(F[i] - 0.5*v[i]*(Pi/V),2.0);
		fcheck = fcheckFrc + pow(P - P0,2.0);
		fcheckFrc = sqrt(fcheckFrc/vertDOF);
		fcheck = sqrt(fcheck/vertDOF);

		// update iterator
		fireit++;
	}
	// check if FIRE converged
	if (fireit == itmax){
		cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
		exit(1);
	}
	else{
		cout << endl;
		cout << "===========================================" << endl;
		cout << " 	F I R E 								" << endl;
		cout << " 	E N T H A L P Y  						" << endl;
		cout << "	M I N I M I Z A T I O N 				" << endl;
		cout << "	C O N V E R G E D! 						" << endl;
		cout << "===========================================" << endl;
		cout << endl;
		cout << "	** fireit 	= " << fireit << endl;
		cout << "	** fcheck 	= " << fcheck << endl;
		cout << " 	** fcheckFrc 	= " << fcheckFrc << endl;
		cout << "	** U 		= " << U << endl;
		cout << "	** dt 		= " << dt << endl;
		cout << "	** PFIRE 	= " << PFIRE << endl;
		cout << "	** alpha 	= " << alpha << endl;
		cout << "	** npPos 	= " << npPos << endl;
		cout << "	** npNeg 	= " << npNeg << endl;
		cout << "	** V  		= " << V << endl;
		cout << "	** P 		= " << Pinst << endl;
		cout << "	** Pi 		= " << Pi << endl;
		cout << endl << endl;
	}
}







/******************************

	P R I N T I N G

*******************************/


void adcm2D::printADCM2DConfiguration(){
	// local variables
	int ci, cj, vi, gi, ctmp, zc, zv;
	double xi, yi, dx, dy, Lx, Ly;

	// check if pos object is open
	if (!posout.is_open()) {
		cerr << "** ERROR: in printADCM2DConfiguration, posout is not open, but function call will try to use. Ending here." << endl;
		exit(1);
	}
	else
		cout << "** In printADCM2DConfiguration, printing particle positions to file..." << endl;

	// save box sizes
	Lx = L.at(0);
	Ly = L.at(1);

	// print information starting information
	posout << setw(w) << left << "NEWFR"
		   << " " << endl;
	posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << endl;
	posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction2D() << endl;

	// print box sizes
	posout << setw(w) << left << "BOXSZ";
	posout << setw(wnum) << setprecision(pnum) << left << Lx;
	posout << setw(wnum) << setprecision(pnum) << left << Ly;
	posout << endl;

	// print coordinate for rest of the cells
	for (ci = 0; ci < NCELLS; ci++) {
		// cell information
		posout << setw(w) << left << "CINFO";
		posout << setw(w) << left << nv.at(ci);
		posout << setw(wnum) << left << a0.at(ci);
		posout << setw(wnum) << left << area(ci);
		posout << setw(wnum) << left << perimeter(ci);
		posout << endl;

		// get initial vertex positions
		gi = gindex(ci, 0);
		xi = x.at(NDIM * gi);
		yi = x.at(NDIM * gi + 1);

		posout << setw(w) << left << "VINFO";
		posout << setw(w) << left << ci;
		posout << setw(w) << left << 0;

		// output initial vertex information
		posout << setw(wnum) << setprecision(pnum) << right << xi;
		posout << setw(wnum) << setprecision(pnum) << right << yi;
		posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << st.at(gi);
		posout << endl;

		// vertex information for next vertices
		for (vi = 1; vi < nv.at(ci); vi++) {
			// get global vertex index for next vertex
			gi++;

			// get next vertex positions
			dx = x.at(NDIM * gi) - xi;
			if (pbc[0])
				dx -= Lx * round(dx / Lx);
			xi += dx;

			dy = x.at(NDIM * gi + 1) - yi;
			if (pbc[1])
				dy -= Ly * round(dy / Ly);
			yi += dy;

			// Print indexing information
			posout << setw(w) << left << "VINFO";
			posout << setw(w) << left << ci;
			posout << setw(w) << left << vi;

			// output vertex information
			posout << setw(wnum) << setprecision(pnum) << right << xi;
			posout << setw(wnum) << setprecision(pnum) << right << yi;
			posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
			posout << setw(wnum) << setprecision(pnum) << right << st.at(gi);
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR"
		   << " " << endl;
}