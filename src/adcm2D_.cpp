/*

	MAIN FUNCTION DEFINITIONS for ADCM2D class

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
adcm2D::adcm2D(int numcells, int numverts, double sizeDisp, double phi0, double boxLengthScale, double clScale, double initGamma0, int seed) : dpm(numcells, 2, seed) {
	// initialize variables for adcm2D
	useRepulsiveForce();
	useIndependentShapeForce();

	// initialize particles with gaussian sizes
	gaussian2D(sizeDisp, 1.0, numverts);
	st.resize(NVTOT);
	dp.resize(NVTOT);
	fill(st.begin(), st.end(), initGamma0);
	fill(dp.begin(), dp.end(), 0.0);
	gam0 = initGamma0;
	Posm = 0.0;
	W = 0.0;

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

	// initialize surface tension matrix 
	stMat.resize(NCELLS+1);
	for (int ci=0; ci<NCELLS+1; ci++){
		stMat.at(ci).resize(NCELLS+1);
		fill(stMat.at(ci).begin(), stMat.at(ci).end(), gam0);
	}
}


// constructor for reading in configuration from file
adcm2D::adcm2D(string &inputFileStr, int seed, double boxLengthScale) : dpm(2) {
	// open file
	ifstream inputobj(inputFileStr.c_str());
	if (!inputobj.is_open()){
		cerr << "** ERROR: In meso2D constructor, could not open file " << inputFileStr << ", ending here. " << endl;
		exit(1);
	}

	// initialize variables for adcm2D
	useRepulsiveForce();
	useIndependentShapeForce();

	// local variables
	int nvtmp, ci, vi, i;
	double val;
	double lxtmp, lytmp;
	double a0tmp, gamtmp;
	double xtmp, ytmp, rtmp;
	string inputStr;

	// LINE 1: should be NEWFR
	getline(inputobj, inputStr);
	if (inputStr.compare(0,5,"NEWFR") != 0){
		cerr << "** ERROR: In meso2D constructor, first line of input file NOT NEWFR, first line = " << inputStr << ". Ending." << endl;
		exit(1);
	}

	// read in simulation information from header
	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"NUMCL %d",&NCELLS);
	cout << "\t ** " << inputStr << endl;

	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"PACKF %lf",&val);
	cout << "\t ** " << inputStr << endl;

	// initialize box lengths
	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"BOXSZ %lf %lf",&lxtmp,&lytmp);
	cout << "\t ** " << inputStr << endl;

	L.at(0) = lxtmp;
	L.at(1) = lytmp;

	// szList and nv (keep track of global vertex indices)
	nv.resize(NCELLS);
	szList.resize(NCELLS);
	a0.resize(NCELLS);

	// initialize NVTOT to 0
	NVTOT = 0;

	// loop over cells, read in coordinates
	cout << "\n\n** LOOPING OVER DATA, PRINTING INFO..." << endl;
	gam0 = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		// first parse cell info
		getline(inputobj, inputStr);
		sscanf(inputStr.c_str(),"CINFO %d %lf %*lf %*lf",&nvtmp,&a0tmp);
		nv.at(ci) = nvtmp;
		a0.at(ci) = a0tmp;		

		// increment szList, NVTOT
		NVTOT += nvtmp;
		if (ci > 0)
			szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

		// print to console
		cout << "\t ** " << inputStr << endl;

		// loop over vertices, store coordinates
		for (vi=0; vi<nvtmp; vi++){
			// parse vertex coordinate info
			getline(inputobj, inputStr);
			sscanf(inputStr.c_str(),"VINFO %*d %*d %lf %lf %lf %lf",&xtmp,&ytmp,&rtmp,&gamtmp);

			// print to console
			cout << "\t ** " << inputStr << endl;

			// push back 
			x.push_back(xtmp);
			x.push_back(ytmp);
			r.push_back(rtmp);
			st.push_back(gamtmp);
			gam0 += gamtmp;
		}
	}
	vertDOF = NDIM * NVTOT;
	gam0 /= vertDOF;
	targetLength = (2.0 * sqrt(PI)) / (NVTOT / NCELLS) ;
	cout << "** NVTOT = " << NVTOT << ", vertDOF = " << vertDOF << ", NCELLS = " << NCELLS << endl;
	cout << "** FINISHED LOOPING OVER DATA, CLOSING INPUT FILE OBJ\n" << endl;

	// initialize vertex indexing
	initializeVertexIndexing2D();

	// other variables that need to be initialized
	setka(1.0);
	setkl(0.0);
	setkb(0.0);
	setkc(0.1);
	setl1(0.0);
	setl2(0.0);

	dp.resize(NVTOT);
	fill(dp.begin(), dp.end(), 0.0);
	Posm = 0.0;
	W = 0.0;

	// initialize surface tension matrix
	stMat.resize(NCELLS+1);
	for (int ci=0; ci<NCELLS+1; ci++){
		stMat.at(ci).resize(NCELLS+1);
		fill(stMat.at(ci).begin(), stMat.at(ci).end(), 0);
	}

	// initialize neighbor linked list
	initializeNeighborLinkedList2D(boxLengthScale);

	// close input file object
	inputobj.close();

	// seed random number generator
	srand48(seed);
}


// set all a0 to be regular values for osmotic pressure
void adcm2D::regularizeA0(){
	// local variables
	int ci;
	double acurr, dtheta, trigfactor, lbar; 

	// loop over cells, set a0
	for (ci=0; ci<NCELLS; ci++){
		// get current area
		acurr = area(ci);

		// angle
		dtheta = (2.0 * PI) / nv.at(ci);

		// trig factor
		trigfactor = (1.0 - cos(dtheta)) / sin(dtheta);

		// segment length
		lbar = perimeter(ci) / nv.at(ci);

		// set a0
		a0.at(ci) = acurr - (lbar / ((Posm * lbar) + 2.0 * gam0 * trigfactor));
	}
}


/******************************

	C I R C U L O - L I N E

	G E O M E T R Y

*******************************/


// get cell "center-of-segment"
void adcm2D::cellCoS(int ci, double &cx, double &cy){
	// local variables
	int vi, gi = szList.at(ci);
	double xcurr, ycurr, dx, dy, l;

	// stay in box with vertex 0
	xcurr = x[NDIM * gi];
	ycurr = x[NDIM * gi + 1];

	// loop over vertices
	cx = 0.0;
	cy = 0.0;
	for (vi=0; vi<nv[ci]; vi++){
		// get distance to next segment 
		dx = x[NDIM * ip1[gi]] - xcurr;
		dy = x[NDIM * ip1[gi] + 1] - ycurr;
		if (pbc[0])
			dx -= L[0] * round(dx / L[0]);
		if (pbc[1])
			dy -= L[1] * round(dy / L[1]);

		// seg length
		l = sqrt(dx * dx + dy * dy);
		
		// weight cx, cy
		cx += xcurr * l;
		cy += ycurr * l;

		// update current position
		xcurr += dx;
		ycurr += dy;

		// increment gi
		gi++;
	}

	// scale by perimeter
	cx /= perimeter(ci);
	cy /= perimeter(ci);
}

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

// points from edge to vertex (PASS EDGE IN FIRST)
double adcm2D::edge2VertexDistance(const int ge, const int gv, double &hx, double &hy, double &tev){
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
				(this->*pwFrc)(gi, gj);

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

					// cell index of gj
					cindices(cj, vj, gj);

					if (ci == cj) {
						pj = list[pj];
						continue;
					}

					// compute force between circulolines
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
void adcm2D::SRRepulsivePWForce(const int gi, const int gj, bool& vivj, bool &viej, bool &vjei){
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
	if (dlijApprox < 4.0 * lij) {
		// projection of vertex j onto edge i
		hr_j2i = edge2VertexDistance(gi, gj, hx_j2i, hy_j2i, tij);

		// projection of vertex i onto edge j
		hr_i2j = edge2VertexDistance(gj, gi, hx_i2j, hy_i2j, tji);
		
		// check determine force based on projection
		if (tij < 0 && dr < sij){
			// report 
			vivj = 1;

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

			// update surface tension
			st.at(gi) = gam0 * (1.0 - W);
			st.at(gj) = gam0 * (1.0 - W);
		}
		else if (tij > 0 && tij < 1 && hr_j2i < sij) {
			// report
			vjei = 1;

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

			// update surface tension
			st.at(gi) = gam0 * (1.0 - W);
			st.at(gj) = gam0 * (1.0 - W);
		}

		// also check projection from edge j to vertex i (same as tij case, but note use of -dx and -dy)
		if (tji > 0 && tji < 1 && hr_i2j < sij) {
			// report
			viej = 1;

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

			// update surface tension
			st.at(gi) = gam0 * (1.0 - W);
			st.at(gj) = gam0 * (1.0 - W);
		}
	}
}

// short-range (SR) _attractive pairwise force
void adcm2D::SRAttractivePWForce(const int gi, const int gj, bool &vivj, bool &viej, bool &vjei){
	// local variables
	int ci, vi, cj, vj;
	double tij, tji;
	double hr_i2j, hr_j2i, hx_i2j, hx_j2i, hy_i2j, hy_j2i;
	double ftmp, dfx, dfy;
	double xij;

	// get ci and cj (for surface tension update)
	cindices(ci, vi, gi);
	cindices(cj, vj, gj);

	// contact distance
	double sij = r[gi] + r[gj];

	// shell, cutoff distance
	double shellij = (1.0 + l2)*sij;
	double cutij = (1.0 + l1)*sij; 
	double kint = (kc * l1) / (l2 - l1);

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
	if (dlijApprox < 4.0 * lij) {
		// projection of vertex j onto edge i
		hr_j2i = edge2VertexDistance(gi, gj, hx_j2i, hy_j2i, tij);

		// projection of vertex i onto edge j
		hr_i2j = edge2VertexDistance(gj, gi, hx_i2j, hy_i2j, tji);
		
		// check determine force based on projection
		// if (tij < 0 && dr < shellij){
		if (tij < 0 && dr < sij){
			// report 
			vivj = 1;
			
			// force
			xij = dr / sij;
			if (dr > cutij){
				ftmp = kint * (1.0 + l2 - xij) / sij;
				addU(-0.5 * kint * pow(1.0 + l2 - xij, 2.0));
			}
			else{
				ftmp = -kc * (1 - xij) / sij;
				addU(0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2));
			}

			// force elements
			dfx = ftmp * (dx / dr);
			dfy = ftmp * (dy / dr);

			// add to force
			addF(dfx, gi, 0);
			addF(dfy, gi, 1);
			addF(-dfx, gj, 0);
			addF(-dfy, gj, 1);

			// add to pressure
			Pinst += ftmp * (dr / L[0]);
		}
		else if (tij > 0 && tij < 1 && hr_j2i < shellij) {
			// report
			vjei = 1;

			// force info
			xij = hr_j2i / sij;
			if (hr_j2i > cutij){
				ftmp = kint * (1.0 + l2 - xij) / sij;
				addU(-0.5 * kint * pow(1.0 + l2 - xij, 2.0));
			}
			else{
				ftmp = -kc * (1 - xij) / sij;
				addU(0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2));
			}

			// force elements
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
			addU(0.5 * kint * pow(1 - (hr_j2i / sij), 2.0));

			// add to pressure
			Pinst += ftmp * (hr_j2i / L[0]);
		}

		// also check projection from edge j to vertex i (same as tij case, but note use of -dx and -dy)
		if (tji > 0 && tji < 1 && hr_i2j < shellij) {
			// report
			viej = 1;

			// force info
			xij = hr_i2j / sij;
			if (hr_i2j > cutij){
				ftmp = kint * (1.0 + l2 - xij) / sij;
				addU(-0.5 * kint * pow(1.0 + l2 - xij, 2.0));
			}
			else{
				ftmp = -kc * (1 - xij) / sij;
				addU(0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2));
			}

			// force elements
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
			addU(0.5 * kint * pow(1 - (hr_i2j / sij), 2.0));

			// add to pressure
			Pinst += ftmp * (hr_i2j / L[0]);
		}
	}
}


// function to add to total force with adhesive force (vertex-vertex contacts)
// NOTE: vector (dx, dy) points from gv1 to gv2
double adcm2D::vvSoftAdhesionForce(const int gv1, const int gv2, const double dr, const double dx, const double dy){
	double dfx = 0.0, dfy = 0.0;
	return vvSoftAdhesionForce(gv1, gv2, dr, dx, dy, dfx, dfy);
}

double adcm2D::vvSoftAdhesionForce(const int gv1, const int gv2, const double dr, const double dx, const double dy, double &dfx, double &dfy){
	// local variables
	double val, sij, xij, shellij, cutij;
	double ftmp = 0.0;
	double kint = (kc * l1) / (l2 - l1);

	// contact distances
	sij = r[gv1] + r[gv2];
	cutij = (1.0 + l1) * sij;
	shellij = (1.0 + l2) * sij;

	// compute and add to forces, energy, pressure
	xij = dr / sij;
	if (dr > shellij)
		ftmp = 0.0;
	else if (dr < shellij && dr > cutij) {
		ftmp = kint * (1.0 + l2 - xij) / sij;
		addU(-0.5 * kint * pow(1.0 + l2 - xij, 2.0));
	}
	else {
		ftmp = -kc * (1 - xij) / sij;
		addU(0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2));
	}

	// force elements
	dfx = ftmp * (dx / dr);
	dfy = ftmp * (dy / dr);

	// add to force
	addF(dfx, gv1, 0);
	addF(dfy, gv1, 1);
	addF(-dfx, gv2, 0);
	addF(-dfy, gv2, 1);

	// add to pressure
	Pinst += ftmp * (dr / L[0]);

	// return ftmp
	return ftmp;
}


// dx, dy should point from edge to vertex 
double adcm2D::evSoftAdhesionForce(const int ge, const int gv, const double dr, const double dx, const double dy, const double t){
	double dfx = 0.0, dfy = 0.0;
	return evSoftAdhesionForce(ge, gv, dr, dx, dy, t, dfx, dfy);
}

double adcm2D::evSoftAdhesionForce(const int ge, const int gv, const double dr, const double dx, const double dy, const double t, double &dfx, double &dfy){
	// local variables
	double val, sij, xij, shellij, cutij;
	double ftmp = 0.0;
	double kint = (kc * l1) / (l2 - l1);

	// contact distances
	sij = r[ge] + r[gv];
	cutij = (1.0 + l1) * sij;
	shellij = (1.0 + l2) * sij;

	// compute and add to forces, energy, pressure
	xij = dr / sij;
	if (dr > shellij)
		ftmp = 0.0;
	else if (dr < shellij && dr > cutij) {
		ftmp = kint * (1.0 + l2 - xij) / sij;
		addU(-0.5 * kint * pow(1.0 + l2 - xij, 2.0));
	}
	else {
		ftmp = -kc * (1 - xij) / sij;
		addU(0.5 * kc * (pow(1.0 - xij, 2.0) - l1 * l2));
	}

	// force elements
	dfx = ftmp * (dx / dr);
	dfy = ftmp * (dy / dr);

	// add to force
	addF(dfx * (1 - t), ge, 0);
	addF(dfy * (1 - t), ge, 1);
	addF(dfx * t, ip1[ge], 0);
	addF(dfy * t, ip1[ge], 1);
	addF(-dfx, gv, 0);
	addF(-dfy, gv, 1);

	// add to pressure
	Pinst += ftmp * (dr / L[0]);

	// return ftmp
	return ftmp;
}


// shape forces for systems with only surface tension (ASSUMING kb, kl == 0)
// TO-DO: NEED TO UPDATE AREA FORCE TO REFLECT P ~ 1/(A - A0) relation
// Look for change to factor of 2 in tension
void adcm2D::adcm2DShapeForces(){
	// local variables
	int ci, gi, vi, nvtmp;
	double fa, cx, cy, xi, yi;
	double a0tmp, atmp, sti, stim1;
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
				stim1 = max(st[im1[gi]], 0.0);

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

		// segment info
		lix = rip1x - rix;
		liy = rip1y - riy;
		li = sqrt(lix * lix + liy * liy);
		sti = max(st[gi], 0.0);

		// add to force
		F[NDIM * gi] += (sti * (lix / li)) - (stim1 * (lim1x / lim1));
		F[NDIM * gi + 1] += (sti * (liy / li)) - (stim1 * (lim1y / lim1));

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
		stim1 = sti;
	}
}

// shape forces for systems with surface tension & contractility (l0 = l0reg - gam/kl)
void adcm2D::adcm2DContractileForces(){
	// local variables
	int ci, gi, vi, nvtmp;
	double fa, cx, cy, xi, yi;
	double a0tmp, atmp, l0i, l0im1, l0reg;
	double dx, dy, da, dli, dlim1;
	double lim1x, lim1y, lix, liy, lip1x, lip1y, li = 0.0, lim1 = 0.0;
	double rim1x, rim1y, rix, riy, rip1x, rip1y;

	// loop over vertices, add to force
	ci = 0;
	for (gi = 0; gi < NVTOT; gi++) {
		// -- Area force (and get cell index ci)
		if (ci < NCELLS) {
			if (gi == szList.at(ci)) {
				// shape information
				nvtmp = nv.at(ci);
				a0tmp = a0.at(ci);

				// effect preferred segment length
				l0reg = sqrt(4.0 * PI * a0tmp) / nvtmp;

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

		// -- Contractility force

		// segment lengths
		lix 	= rip1x - rix;
		liy		= rip1y - riy;
		li 		= sqrt(lix * lix + liy * liy);

		// segment strains from preferred lengths
		l0im1 	= l0reg;
		l0i 	= l0reg;

		dlim1  	= (lim1/l0im1) - 1.0;
		dli 	= (li/l0i) - 1.0;

		// add to force
		// cout << "gi=" << gi << ",   dli=" << dli << ",    dlim1=" << dlim1 << ",    flx (im1 comp.)=" << -(dlim1 * (lim1x / (lim1 * l0im1))) << ",    fly (im1 comp.)=" << -(dlim1 * (lim1y / (lim1 * l0im1))) << ",   flx=" << (dli * (lix / (li * l0i))) - (dlim1 * (lim1x / (lim1 * l0im1))) << ",   fly=" << (dli * (liy / (li * l0i))) - (dlim1 * (lim1y / (lim1 * l0im1))) << endl;
		F[NDIM * gi] += dli * (lix / (li * l0i)) - dlim1 * (lim1x / (lim1 * l0im1));
		F[NDIM * gi + 1] += dli * (liy / (li * l0i)) - dlim1 * (lim1y / (lim1 * l0im1));

		// add to potential energy
		U += 0.5 * dli * dli;

		// add to pressure and shear stress (TO DO: add shear stress)
		// NOTE: this is not yet the pressure, just the contribution from \partial U / \partial L
		Pinst += (dli * li)/(l0i * L[0]);

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


// shape forces with osmotic pressure from external environment, other cells
// NOTE: NEED TO RECOMPUTE CONTRIBUTION TO EXTERNAL PRESSURE FOR NEW PRESSURE FORM
void adcm2D::adcm2DOsmoticPressureShapeForces(){
	// local variables
	int ci, gi, vi, nvtmp;
	double fa, cx, cy, xi, yi;
	double a0tmp, atmp, sti, stim1, dpi, dpim1;
	double dx, dy, da, dli, dlim1;
	double lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
	double rim1x, rim1y, rix, riy, rip1x, rip1y;

	// set all pressure differences to be with void
	gi = 0;
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<nv.at(ci); vi++){
			// set dp for vertex gi based on area of ci
			dp.at(gi) = Posm - (ka / (area(ci) - a0.at(ci)));
			gi++;
		}
	}

	// loop over vertices, add to force
	ci = 0;
	for (gi = 0; gi<NVTOT; gi++) {
		// -- Area force (and get cell index ci)
		if (ci < NCELLS) {
			if (gi == szList[ci]) {
				// shape information
				nvtmp = nv.at(ci);

				// update potential energy, pressure
				da 		= area(ci) - a0.at(ci);
				U 		= Posm * area(ci) - ka * log(da);
				Pinst 	+= ((2.0 * area(ci))/L[0]) * (Posm - (ka/da));

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
				stim1 = max(st[im1[gi]], 0.0);
				dpim1 = dp[im1[gi]];

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


		// segment & pressure info
		lix = rip1x - rix;
		liy = rip1y - riy;
		li = sqrt(lix * lix + liy * liy);
		sti = max(st[gi], 0.0);
		dpi = dp[gi];

		// -- Area force
		F[NDIM * gi] += 0.5 * (dpim1 * rim1y - dpi * rip1y);
		F[NDIM * gi + 1] += 0.5 * (dpi * rip1x - dpim1 * rim1x);

		// -- Surface tension force
		F[NDIM * gi] += (sti * (lix / li)) - (stim1 * (lim1x / lim1));
		F[NDIM * gi + 1] += (sti * (liy / li)) - (stim1 * (lim1y / lim1));

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
		stim1 = sti;
		dpim1 = dpi;
	}

	// // reset all pressure differences to be with void (updated in interaction force)
	// gi = 0;
	// for (ci=0; ci<NCELLS; ci++){
	// 	for (vi=0; vi<nv[ci]; vi++){
	// 		// set dp for vertex gi based on area of ci
	// 		dp.at(gi) = Posm - (ka / (area(ci) - a0[ci]));
	// 		gi++;
	// 	}
	// }
}


/******************************

	B O N D  &  V E R T E X

	D Y N A M I C S

*******************************/

// check vertices to see if some should be deleted or added
// TO DO: 
// * Change to removing via segments, not vertices
// * Add single criterion, to either add or remove based on density consideration and not single segment size
// * OR change criteria to be based on local curvature ... min length changes based on curvature radius (more curvature, need more vertices)
void adcm2D::checkVertices(){
	// local variables
	int gi;
	double li;

	// frac for addition / deletion
	const double changeFrac = 0.2;
	const double addFrac = 2.0;
	const double delFrac = 0.5; 

	// lists for deletion and addition
	vector<int> verts2Del;
	vector<int> verts2Add;

	// loop over all vertices, check length, delete first
	for (gi=0; gi<NVTOT; gi++){
		// get segment length
		li = seg(gi);
		if (li < delFrac * targetLength)
			verts2Del.push_back(gi);
	}

	// loop over vertices to delete, delete them
	for (gi=0; gi<verts2Del.size(); gi++)
		deleteVertex(verts2Del.at(gi) - gi);

	// loop over all vertices 
	for (gi=0; gi<NVTOT; gi++){
		// get new segment length
		li = seg(gi);
		if (li > addFrac * targetLength)
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


// Energy minimization with a central potential centered at origin
void adcm2D::ctrPotentialFIRE(double Ftol, double kctr, double dt0, bool printRelaxation){
	// local variables
	int i, d, gi;
	double cx, cy, dcx, dcy, dfx, dfy, ctrRMSD, cx0, cy0;
	double shpTrq, shpF, intTrq, intF, ctrTrq, ctrF;

	// check to see if cell linked-list has been initialized
	if (NBX == -1){
		cerr << "	** ERROR: In adcm2D::mesoEnthalpyFIRE, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
		exit(1);
	}

	// FIRE variables
	double PFIRE, fnorm, fcheck, dPcheck, vnorm, alpha, dtmax, dtmin;
	int npPos, npNeg, fireit;

	// set dt
	dt = dt0;

	// Initialize FIRE variables
	PFIRE  		= 0;	
	fnorm 		= 0;
	vnorm 		= 0;
	alpha   	= alpha0;

	dtmax   	= 2.0 * dt;
	dtmin   	= 0.1 * dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;

	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);
	adcm2DForceUpdate();

	// DEBUGGING: print out forces
	string forcef = "measureFriction.frc";
    ofstream forceOuputObj(forcef.c_str());

	// relax forces using FIRE
	while ((fcheck > Ftol || fireit < NDELAY) && fireit < itmax){
	// while ((fcheck > Ftol || fireit < NDELAY) && fireit < 20000){
		// VV VELOCITY UPDATE #1
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*dt*F[i];

		// compute PFIRE
		PFIRE = 0.0;
		for (i=0; i<vertDOF; i++)
			PFIRE += v[i]*F[i];

		// print to console
		if (fireit % NSKIP == 0){
		// if (fireit > 19000){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << " 	F I R E 								" << endl;
			cout << "	M I N I M I Z A T I O N 				" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit 	= " << fireit << endl;
			cout << "	** fcheck 	= " << fcheck << endl;
			cout << "	** U 		= " << U << endl;
			cout << "	** dt 		= " << dt << endl;
			cout << "	** PFIRE 	= " << PFIRE << endl;
			cout << "	** alpha 	= " << alpha << endl;
			cout << "	** npPos 	= " << npPos << endl;
			cout << "	** npNeg 	= " << npNeg << endl;
			cout << "	** cx 		= " << setprecision(12) << cx << endl;
			cout << "	** cy 		= " << setprecision(12) << cy << endl;
			cout << " 	** dc 		= " << sqrt(pow(0.5 * L[0] - cx, 2.0) + pow(0.5 * L[1] - cy, 2.0)) << endl;

			// print out torque information
			fill(F.begin(), F.end(), 0.0);

			// shape forces
			(*this.*shpFrc)();
			shpTrq = 0.0;
			shpF = 0.0;
			for (gi=0; gi<NVTOT; gi++){
				shpTrq += x[NDIM * gi] * F[NDIM * gi + 1] - x[NDIM * gi + 1] * F[NDIM * gi];
				shpF += F[NDIM * gi] * F[NDIM * gi] + F[NDIM * gi + 1] * F[NDIM * gi + 1];
			}
			shpF = sqrt(shpF / NVTOT);
			cout << "	** shpTrq 	= " << setprecision(12) << shpTrq << endl; 
			cout << "	** shpF 	= " << setprecision(12) << shpF << endl; 

			// interaction forces
			circuloLinePWForceUpdate();
			intTrq = 0.0;
			intF = 0.0;
			for (gi=0; gi<NVTOT; gi++){
				intTrq += x[NDIM * gi] * F[NDIM * gi + 1] - x[NDIM * gi + 1] * F[NDIM * gi];
				intF += F[NDIM * gi] * F[NDIM * gi] + F[NDIM * gi + 1] * F[NDIM * gi + 1];
			}
			intF = sqrt(intF / NVTOT);
			cout << "	** intTrq 	= " << setprecision(12) << shpTrq << endl; 
			cout << "	** intF 	= " << setprecision(12) << shpF << endl;

			// central potential
			for (int ci=0; ci<NCELLS; ci++){
				// get centers of mass
				com2D(ci, cx, cy);

				// compute force
				dcx = 0.5 * L[0] - cx;
				dcy = 0.5 * L[1] - cy;

				// add force to center of mass if far enough away
				dfx = kctr * dcx;
				dfy = kctr * dcy;
					
				addComF(ci, 0, dfx);
				addComF(ci, 1, dfy);
			}
			ctrTrq = 0.0;
			ctrF = 0.0;
			for (gi=0; gi<NVTOT; gi++){
				ctrTrq += x[NDIM * gi] * F[NDIM * gi + 1] - x[NDIM * gi + 1] * F[NDIM * gi];
				ctrF += F[NDIM * gi] * F[NDIM * gi] + F[NDIM * gi + 1] * F[NDIM * gi + 1];
			}
			ctrF = sqrt(ctrF / NVTOT);
			cout << "	** ctrTrq 	= " << setprecision(12) << ctrTrq << endl; 
			cout << "	** ctrF 	= " << setprecision(12) << ctrF << endl; 


			// print it
			if (printRelaxation && posout.is_open()){
				printADCM2DConfiguration();
				printInstantaneousForces(forceOuputObj);
			}
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
				x[i] -= 0.5*dt*v[i];

				// reset vertex velocities
				v[i] = 0.0;
			}

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
			fnorm 	+= pow(F[i],2.0);
			vnorm 	+= v[i]*v[i];
		}
		fnorm = sqrt(fnorm);
		vnorm = sqrt(vnorm);

		// update velocities (s.d. vs inertial dynamics) only if forces are acting
		if (fnorm > 0){
			for (i=0; i<vertDOF; i++)
				v[i] = (1 - alpha)*v[i] + alpha*(F[i]/fnorm)*vnorm;
		}

		// VV POSITION UPDATE
		for (i=0; i<vertDOF; i++)
			x[i] += dt*v[i];

		// update forces, pressure
		adcm2DForceUpdate();

		// add central potential
        for (int ci=0; ci<NCELLS; ci++){
            // get centers of mass
            com2D(ci, cx, cy);

            // compute force
            dcx = 0.5 * L[0] - cx;
            dcy = 0.5 * L[1] - cy;
			// dcx = 0.5 * L[0] - x[NDIM * szList.at(ci)];
			// dcy = 0.5 * L[1] - x[NDIM * szList.at(ci) + 1];

            // add force to center of mass if far enough away
            dfx = kctr * dcx;
            dfy = kctr * dcy;
                
            addComF(ci, 0, dfx);
            addComF(ci, 1, dfy);
			// F[NDIM * szList.at(ci)] += dfx;
			// F[NDIM * szList.at(ci) + 1] += dfy;
        }

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*dt*F[i];

		// update fcheck based on fnorm (= force per degree of freedom)
		fcheck = 0.0;
		for (i=0; i<vertDOF; i++)
			fcheck += pow(F[i],2.0);
		fcheck = sqrt(fcheck/vertDOF);

		// update iterator
		fireit++;
	}

	// DEBUGGING: close force ouput file
	forceOuputObj.close();

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
		cout << "	** U 		= " << U << endl;
		cout << "	** dt 		= " << dt << endl;
		cout << "	** PFIRE 	= " << PFIRE << endl;
		cout << "	** alpha 	= " << alpha << endl;
		cout << "	** npPos 	= " << npPos << endl;
		cout << "	** npNeg 	= " << npNeg << endl;
		cout << "	** P 		= " << Pinst << endl;
		cout << endl << endl;
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
		for (d = 0; d < NDIM; d++){
			sb[d] = floor(L[d] / lb[d]);
			if (sb[d] < 3)
				sb[d] = 3;
		}
			

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


// Enthalpy minimization, different dt for P and degrees of freedom
void adcm2D::nphSplitFIRE(double Ftol, double dPtol, double P0, double dt0, bool printCompression){
	// local variables
	int i, d;

	// DEBUGGING: print out forces
	string forcef = "measureFriction.frc";
    ofstream forceOuputObj(forcef.c_str());

	// check whether neighbor list changes
	bool boxNumChange = false;

	// box size, momentum, and internal virial pressure
	Pinst = 0.0;
	double V=L[0]*L[1], Pi=0.0, P=0.0;
	double pdt, bdt;

	// FIRE variables for particles
	double pP, pFNRM, pVNRM, pAL, FCHECK;
	int pNPPOS, pNPNEG;

	// FIRE variables for box
	double bdtmax, bdtmin, bP, bFNRM, bVNRM, bAL, DPCHECK;
	int bNPPOS, bNPNEG;

	// FIRE variables for everything
	double dtmax, dtmin;
	int fireit;

	// set dt based on geometric parameters
	pdt = dt0;
	bdt = dt0;

	// Initialize FIRE variables
	pAL   		= alpha0;
	bAL 		= alpha0;

	dtmax   	= 5.0 * dt0;
	dtmin   	= 0.2 * dt0;

	bdtmax 		= 5.0 * bdt;
	bdtmin 		= 0.2 * bdt;

	pNPPOS    	= 0;
	pNPNEG      = 0;
	bNPPOS 		= 0;
	bNPNEG 	 	= 0;

	fireit    	= 0;
	FCHECK  	= 10*Ftol;
	DPCHECK 	= 10*dPtol;

	// initialize loop
	while ((FCHECK > Ftol || DPCHECK > dPtol || fireit < NDELAY) && fireit < itmax){
		// VV VELOCITY UPDATE #1
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5 * pdt * (F[i] - 0.5 * v[i] * (Pi / V));
		Pi += 0.5 * bdt * (P - P0);

		// compute projection of force onto velocity
		pP = 0.0;
		for (i=0; i<vertDOF; i++)
			pP += v[i] * (F[i] - 0.5 * v[i] * (Pi / V));
		bP = Pi * (P - P0);

		// print to console
		if (fireit % 5000 == 0){
			cout << endl << endl;
			cout << "===============================" << endl;
			cout << " 	S P L I T  F I R E 			" << endl;
			cout << " 	E N T H A L P Y  			" << endl;
			cout << "	M I N I M I Z A T I O N 	" << endl;
			cout << "================================" << endl;
			cout << endl;
			cout << "	** fireit 	= " << fireit << endl;
			cout << "	** FCHECK 	= " << FCHECK << endl;
			cout << " 	** DPCHECK 	= " << DPCHECK << endl;
			cout << "	** pdt 		= " << pdt << endl;
			cout << "	** bdt 		= " << bdt << endl;
			cout << "	** V  		= " << V << endl;
			cout << "	** P 		= " << P << endl;
			cout << "	** Pi 		= " << Pi << endl << endl;
			
			cout << "	Particles:" << endl;
			cout << "	** pP 		= " << pP << endl;
			cout << "	** pAL 		= " << pAL << endl;
			cout << "	** pFNRM 	= " << pFNRM << endl;
			cout << "	** pVNRM 	= " << pVNRM << endl;
			cout << "	** pNPPOS 	= " << pNPPOS << endl;
			cout << "	** pNPNEG 	= " << pNPNEG << endl;

			cout << "	Box:" << endl;
			cout << "	** bP 		= " << bP << endl;
			cout << "	** bAL 		= " << bAL << endl;
			cout << "	** bFNRM 	= " << bFNRM << endl;
			cout << "	** bVNRM 	= " << bVNRM << endl;
			cout << "	** bNPPOS 	= " << bNPPOS << endl;
			cout << "	** bNPNEG 	= " << bNPNEG << endl;
			
			// print it
			if (printCompression && posout.is_open()){
				printADCM2DConfiguration();
				// printInstantaneousForces(forceOuputObj);
			}
		}


		// -- ADJUST PARTICLE DYNAMICS

		// Adjust simulation based on net motion of degrees of freedom
		if (pP > 0){
			// counters
			pNPPOS++;
			pNPNEG = 0;

			// alter simulation if enough positive steps have been taken
			if (pNPPOS > NDELAY){
				// change time step
				if (pdt*finc < dtmax)
					pdt *= finc;

				// decrease alpha
				pAL *= falpha;
			}
		}
		else{
			// couters
			pNPPOS = 0;
			pNPNEG++;

			// check if simulation is stuck
			if (pNPNEG > NNEGMAX){
				cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
				exit(1);
			}

			// take half step backwards, reset velocities
			for (i=0; i<vertDOF; i++){
				x[i] -= 0.5 * pdt * (v[i] + 0.5*x[i]*(Pi/V));
				v[i] = 0.0;
			}

			// decrease time step if past initial delay
			if (fireit > NDELAY){
				// decrease time step 
				if (pdt * fdec > dtmin)
					pdt *= fdec;

				// reset alpha
				pAL = alpha0;
			}
		}

		// compute fnorm and vnorm for particles
		pFNRM = 0.0;
		pVNRM = 0.0;
		for (i=0; i<vertDOF; i++){
			pFNRM 	+= F[i] * F[i];
			// pFNRM 	+= pow((F[i] - 0.5 * v[i] * (Pi / V)), 2.0);
			pVNRM 	+= v[i] * v[i];
		}
		pFNRM = sqrt(pFNRM);
		pVNRM = sqrt(pVNRM);


		// adjust velocities
		if (pFNRM > 0){
			for (i=0; i<vertDOF; i++){
				// v[i] = (1 - pAL) * v[i] + pAL * ((F[i] - 0.5 * v[i] * (Pi / V))/ pFNRM) * pVNRM;
				v[i] = (1 - pAL) * v[i] + pAL * (F[i]/ pFNRM) * pVNRM;
			}
				
		}



		// -- ADJUST BOX DYNAMICS

		// Adjust simulation based on net motion of degrees of freedom
		if (bP > 0){
			// counters
			bNPPOS++;
			bNPNEG = 0;

			// alter simulation if enough positive steps have been taken
			if (bNPPOS > NDELAY){
				// change time step
				if (bdt*finc < bdtmax)
					bdt *= finc;

				// decrease alpha
				bAL *= falpha;
			}
		}
		else{
			// couters
			bNPPOS = 0;
			bNPNEG++;

			// check if simulation is stuck
			if (bNPNEG > NNEGMAX){
				cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
				exit(1);
			}

			// take half step backwards, reset velocities
			V -= 0.5 * bdt * Pi;

			// decrease time step if past initial delay
			if (fireit > NDELAY){
				// decrease time step 
				if (bdt * fdec > bdtmin)
					bdt *= fdec;

				// reset alpha
				bAL = alpha0;
			}
		}

		// compute fnorm and vnorm for particles
		bFNRM = abs(P - P0);
		bVNRM = abs(Pi);

		// adjust velocities
		if (bFNRM > 0)
			Pi = (1 - bAL) * Pi + bAL * ((P - P0)/bFNRM) * bVNRM;



		// VV POSITION UPDATE
		for (i=0; i<vertDOF; i++)
			x[i] += pdt * (v[i] + 0.5 * x[i] * (Pi/V));
		V += 10. * bdt * Pi;
		L[0] = sqrt(V);
		L[1] = sqrt(V);

		// update NN based on new box size
		boxNumChange = false;
		for (d = 0; d < NDIM; d++){
			if (round(L[d] / lb[d]) != sb[d]){
				// update number of boxes
				sb[d] = round(L[d] / lb[d]);
				if (sb[d] < 3)
					sb[d] = 3;
					
				// reset box size
				lb[d] = L[d] / sb[d];

				// change bool for reset
				boxNumChange = true;
			}
		}

		// reset nn list if number of boxes change
		if (boxNumChange){
			cerr << "resetting NNs" << endl;
			resetNeighborLinkedListCellNNs();
			printNeighborList();
		}

		// update forces, pressure
		adcm2DForceUpdate();
		P = Pinst;

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5 * pdt * (F[i] - 0.5 * v[i] * (Pi / V));
		Pi += 0.5 * bdt * (P - P0);
		// Pi += bdt * (P - P0); // dropping 0.5 means that mass is half as large

		// update fcheck based on fnorm (= force per degree of freedom)
		FCHECK = 0.0;
		for (i=0; i<vertDOF; i++)
			FCHECK += F[i] * F[i];
		FCHECK = sqrt(FCHECK / vertDOF);
		DPCHECK = abs(P - P0);

		// update iterator
		fireit++;
	}
	forceOuputObj.close();
	// check if FIRE converged
	if (fireit == itmax){
		cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
		exit(1);
	}
	else{
		cout << endl << endl;
			cout << "===============================" << endl;
			cout << " 	S P L I T  F I R E 			" << endl;
			cout << " 	E N T H A L P Y  			" << endl;
			cout << "	M I N I M I Z A T I O N 	" << endl;
			cout << "================================" << endl;
			cout << endl;
			cout << "	** fireit 	= " << fireit << endl;
			cout << "	** FCHECK 	= " << FCHECK << endl;
			cout << " 	** DPCHECK 	= " << DPCHECK << endl;
			cout << "	** pdt 		= " << pdt << endl;
			cout << "	** bdt 		= " << bdt << endl;
			cout << "	** V  		= " << V << endl;
			cout << "	** P 		= " << P << endl;
			cout << "	** Pi 		= " << Pi << endl << endl;
			
			cout << "	Particles:" << endl;
			cout << "	** pP 		= " << pP << endl;
			cout << "	** pAL 		= " << pAL << endl;
			cout << "	** pFNRM 	= " << pFNRM << endl;
			cout << "	** pVNRM 	= " << pVNRM << endl;
			cout << "	** pNPPOS 	= " << pNPPOS << endl;
			cout << "	** pNPNEG 	= " << pNPNEG << endl << endl;

			cout << "	Box:" << endl;
			cout << "	** bP 		= " << bP << endl;
			cout << "	** bAL 		= " << bAL << endl;
			cout << "	** bFNRM 	= " << bFNRM << endl;
			cout << "	** bVNRM 	= " << bVNRM << endl;
			cout << "	** bNPPOS 	= " << bNPPOS << endl;
			cout << "	** bNPNEG 	= " << bNPNEG << endl;
			
			// print it
			if (printCompression && posout.is_open())
				printADCM2DConfiguration();
	}
	


}


// Active brownian crawlers
void adcm2D::activeBrownianCrawling(const double Tsim, const double Tprint, const double dt0, const double v0, const double Dr, const double Ds){
	// local variables
	double t = 0.0;
	const double vmin = 0.01 * v0;
	double r1, r2, grv, psitmp, dpsi, cx, cy, dx, dy, dr, ftmp;
	int i, k, ci, vi, gi;
	vector<double> psi(NCELLS, 0.0);

	// set time
	setdt_pure(dt0);

	// keep track of print step
	int nprint = 0;

	// initialize psi in random directions
	for (ci=0; ci<NCELLS; ci++)
		psi.at(ci) = 2.0 * PI * drand48();

	// loop
	k = 0;
	while(t < Tsim){
		// compute forces
		adcm2DForceUpdate();

		// add force based on crawling direction
		gi = 0;
		for (ci=0; ci<NCELLS; ci++){
			// compute forces on all vertices
			com2D(ci, cx, cy);
			for (vi=0; vi<nv[ci]; vi++){
				// relative location to center
				dx = x[NDIM * gi] - cx;
				dx -= L[0] * round(dx / L[0]);

				dy = x[NDIM * gi + 1] - cy;
				dy -= L[1] * round(dy / L[1]);

				dr = sqrt(dx * dx + dy * dy);

				// get psitmp
				psitmp = atan2(dy, dx);

				// get force magnitude
				// dpsi = psitmp - psi.at(ci);
				// dpsi -= 2.0 * PI * round(dpsi / (2.0 * PI));
				// ftmp = v0 * exp(-pow(dpsi, 2.0) / (2.0 * Ds * Ds)) + vmin;

				// add to forces
				// F[NDIM * gi] += ftmp * (dx / dr);
				// F[NDIM * gi + 1] += ftmp * (dy / dr);

				ftmp = (v0 / nv[ci]);
				F[NDIM * gi] += ftmp * cos(psi.at(ci));
				F[NDIM * gi + 1] += ftmp * sin(psi.at(ci));

				// increment gi
				gi++;
			}

			// update psi
			r1 = drand48();
			r2 = drand48();
			grv = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
			psi.at(ci) += sqrt(2.0 * dt * Dr) * grv;
		}

		// Euler update
		for (i=0; i<vertDOF; i++)
			x[i] += dt * F[i];

		// Print
		if (floor(t / Tprint) == nprint){
			cout << endl << endl;
			cout << "===============================" << endl;
			cout << " 	A C T I V E 				" << endl;
			cout << " 	B R O W N I A N 			" << endl;
			cout << " 	C R A W L E R S				" << endl;
			cout << "================================" << endl;
			cout << endl;
			cout << "	** t 		= " << t << endl;
			cout << "	** nprint 	= " << nprint << endl;

			// update nprint
			nprint++;
			
			// print it
			if (posout.is_open())
				printADCM2DConfiguration();
		}

		// update t
		t += dt;
	}
}


// Active tension fluctuations
void adcm2D::activeTensionFluctuations(const double Tsim, const double Tprint, const double dt0, const double deltaST){
// local variables
	double t = 0.0;
	double r1, r2, grv, sttmp, gam0tmp;
	int i, k, ci, cj;

	// set time
	setdt_pure(dt0);

	// noise strength
	const double noiseStrength = deltaST * gam0 * sqrt(dt);

	// keep track of print step
	int nprint = 0;

	// compute mean gam0
	double meanGam0 = 0.0;

	// loop
	k = 0;
	while(t < Tsim){
		// compute forces
		activeTensionForceUpdate();

		// drive fluctuations to surface tensions between cells (or with cells and void)
		meanGam0 = 0.0;
		for (ci=0; ci<NCELLS+1; ci++){
			for (cj=ci; cj<NCELLS+1; cj++){
				// get gam0
				if (ci == 0 || cj == 0)
					gam0tmp = gam0;
				else
					gam0tmp = gam0 * (1.0 - W);
				
				// Gaussian Random Variable using Box-Müller
				r1 = drand48();
				r2 = drand48();
				grv = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
				
				// update surface tension
				sttmp = stMat.at(ci).at(cj);
				// sttmp += dt * (gam0tmp - sttmp) + noiseStrength * grv;
				sttmp += noiseStrength * grv;
				stMat.at(ci).at(cj) = sttmp;
				stMat.at(cj).at(ci) = sttmp;
				meanGam0 += sttmp;
			}
		}
		meanGam0 /= 0.5 * (NCELLS + 2) * (NCELLS + 1);

		// Euler update
		for (i=0; i<vertDOF; i++)
			x[i] += dt * F[i];

		// Print
		if (floor(t / Tprint) == nprint){
			cout << endl << endl;
			cout << "===============================" << endl;
			cout << " 	A C T I V E 				" << endl;
			cout << " 	T E N S I O N 				" << endl;
			cout << " 	F L U C T U A T I O N S		" << endl;
			cout << "===============================" << endl;
			cout << endl;
			cout << "	** t 		= " << t << " / " << Tsim << endl;
			cout << "	** nprint 	= " << nprint << endl;
			cout << "	** meanGam0 = " << meanGam0 << endl;

			// update nprint
			nprint++;
			
			// print it
			if (posout.is_open())
				printADCM2DConfiguration();
		}

		// update t
		t += dt;
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

	// print stresses
	posout << setw(w) << left << "STRSS";
	posout << setw(wnum) << setprecision(pnum) << left << Pinst;
	posout << setw(wnum) << setprecision(pnum) << left << Sinst;
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



// function to print out instantaneous force vectors
void adcm2D::printInstantaneousForces(ofstream &outputObj){
	// first, check if object is open
	if (!outputObj.is_open()) {
		cout << "ERROR: In printInstantaneousForces, outputObj is not open, so ending here." << endl;
		exit(1);
	}
	else {
		cout << "** In printInstantaneousForces, printing forces to file..." << endl;
	}

	// print to start frame
	outputObj << "NEWFR" << endl;

	// local variables
	int ci, vi, cj, vj, gi, gj;
	int gip1, gim1, gjp1, gjm1;
	double sij, shellij, cutij, drijp1, drjip1, lij, dlijApprox;
	double dx, dy, dr, dfx, dfy;
	double ftmpi, ftmpj, trq, tqi, tqj;
	double pxi, pyi, delxi, delyi, tiprint, fxi, fyi;
	double pxj, pyj, delxj, delyj, tjprint, fxj, fyj;

	// projection components from EDGES ((im1, i, ip1), mu) TO VERTEX (j, nu)
	double tim1j, tij, tip1j;
	double hr_ij, hx_ij, hy_ij;
	double hr_im1j, hx_im1j, hy_im1j;
	double hr_ip1j, hx_ip1j, hy_ip1j;

	// projection components from EDGES ((jm1, j, jp1), nu) TO VERTEX (i, mu)
	double tjm1i, tji, tjp1i;
	double hr_ji, hx_ji, hy_ji;
	double hr_jm1i, hx_jm1i, hy_jm1i;
	double hr_jp1i, hx_jp1i, hy_jp1i;

	// for min distance checking
	int pi, pj;
	double dij, dijx, dijy, dji, djix, djiy, ti, tj;
	
	// loop over all pairs of segments, compute force if interacting
	for (gi=0; gi<NVTOT; gi++){
		cindices(ci, vi, gi);
		for (gj=(gi+1); gj<NVTOT; gj++){
			// skip if gj is on same cell as gi
			cindices(cj, vj, gj);
			if (ci == cj)
				continue;

			// reset force counters
			ftmpi = 0.0;
			ftmpj = 0.0;
			
			// -- Compute pairwise force
			
			// contact distances
			sij = r[gi] + r[gj];
			shellij = (1.0 + l2)*sij;
			cutij = (1.0 + l1)*sij; 

			// indices
			gip1 = ip1[gi];
			gim1 = im1[gi];
			gjp1 = ip1[gj];
			gjm1 = im1[gj];

			// distance between vertex i and j
			dx = deltaX(gi, gj, 0);
			dy = deltaX(gi, gj, 1);
			dr = sqrt(dx * dx + dy * dy);

			// distance between vertex i and jp1
			drjip1 = deltaR(gj, gip1);
			drijp1 = deltaR(gi, gjp1);

			// Gut check distance
			lij = 0.5 * (seg(gi) + seg(gj));
		 	dlijApprox = 0.5 * (dr + deltaR(gip1, gjp1));
			if (dlijApprox < 20.0 * lij) {
				// vectors from edge to vertex (EDGE IS PASSED FIRST)
				hr_ij = edge2VertexDistance(gi, gj, hx_ij, hy_ij, tij);
				hr_ji = edge2VertexDistance(gj, gi, hx_ji, hy_ji, tji);
				
				// -- 	force due to projection from edge i onto vertex j, or due to vertex i - vertex j overlap, 
				// 		or due to projection from edge j onto vertex i
				if (dr < shellij) {
					// compute projection from edge im1 to vertex j (PASS EDGE INDEX FIRST)
					hr_im1j = edge2VertexDistance(gim1, gj, hx_im1j, hy_im1j, tim1j);
					hr_jm1i = edge2VertexDistance(gjm1, gi, hx_jm1i, hy_jm1i, tjm1i);
					
					// check projection component from edges on mu to vertex (j, nu), determine relevant minimum distance
					if (tij < 0 && tim1j > 1){
						// vertex - vertex contact
						pi = gi;
						dij = dr;
						dijx = dx;
						dijy = dy;
						ti = -1;
					}
					else if ((tij > 0 && tij < 1) && tim1j > 1){
						// vertex - edge contact
						pi = gi;
						dij = hr_ij;
						dijx = hx_ij;
						dijy = hy_ij;
						ti = tij;
					}
					else if (tij < 0 && (tim1j > 0 && tim1j < 1)){
						pi = gim1;
						dij = hr_im1j;
						dijx = hx_im1j;
						dijy = hy_im1j;
						ti = tim1j;
					}
					else {
						// projection on both i and im1, so need to check minimum
						if (hr_ij < hr_im1j){
							pi = gi;
							dij = hr_ij;
							dijx = hx_ij;
							dijy = hy_ij;
							ti = tij;
						}
						else{
							pi = gim1;
							dij = hr_im1j;
							dijx = hx_im1j;
							dijy = hy_im1j;
							ti = tim1j;
						}
					}

					// do same, but for projections from vertex (i, mu) to edges on cell nu
					if (tji < 0 && tji > 1){
						// vertex - vertex contact
						pj = gj;
						dji = dr;
						djix = -dx;
						djiy = -dy;
						tj = -1;
					}
					else if ((tji > 0 && tji < 1) && tjm1i > 1){
						// vertex - edge contact
						pj = gj;
						dji = hr_ji;
						djix = hx_ji;
						djiy = hy_ji;
						tj = tji;
					}
					else if (tji < 0 && (tjm1i > 0 && tjm1i < 1)){
						pj = gjm1;
						dji = hr_jm1i;
						djix = hx_jm1i;
						djiy = hy_jm1i;
						tj = tjm1i;
					}
					else {
						// projection on both i and im1, so need to check minimum
						if (hr_ji < hr_jm1i){
							pj = gj;
							dji = hr_ji;
							djix = hx_ji;
							djiy = hy_ji;
							tj = tji;
						}
						else{
							pj = gjm1;
							dji = hr_jm1i;
							djix = hx_jm1i;
							djiy = hy_jm1i;
							tj = tjm1i;
						}
					}

					// use force based on ti, tj, minimum distance
					if (dij < dji){
						if (ti < 0){
							ftmpi = vvSoftAdhesionForce(pi, gj, dij, dijx, dijy, dfx, dfy);
							pxi = x[NDIM * gi];
							pyi = x[NDIM * gi + 1];
							delxi = dijx;
							delyi = dijy;
							tiprint = -100;
							fxi = dfx;
							fyi = dfy;
						}
						else {
							ftmpi = evSoftAdhesionForce(pi, gj, dij, dijx, dijy, ti, dfx, dfy);
							pxi = x[NDIM * gi] + ti * segX(pi, 0);
							pyi = x[NDIM * gi + 1] + ti * segX(pi, 1);
							delxi = dijx;
							delyi = dijy;
							tiprint = ti;
							fxi = dfx;
							fyi = dfy;
						}
					}
					else {
						if (tj < 0){
							ftmpj = vvSoftAdhesionForce(pj, gi, dji, djix, djiy, dfx, dfy);
							pxj = x[NDIM * gj];
							pyj = x[NDIM * gj + 1];
							delxj = djix;
							delyj = djiy;
							tjprint = -100;
							fxj = dfx;
							fyj = dfy;
						}
						else {
							ftmpj = evSoftAdhesionForce(pj, gi, dji, djix, djiy, tj, dfx, dfy);
							pxj = x[NDIM * gj] + tj * segX(pj, 0);
							pyj = x[NDIM * gj + 1] + tj * segX(pj, 1);
							delxj = djix;
							delyj = djiy;
							tjprint = tj;
							fxj = dfx;
							fyj = dfy;
						}
					}
				}
				else if (drjip1 > shellij && hr_ij < shellij && tij > 0 && tij < 1) {
					// compute projection onto adjacent edges
					hr_im1j = edge2VertexDistance(gim1, gj, hx_im1j, hy_im1j, tim1j);
					hr_ip1j = edge2VertexDistance(gip1, gj, hx_ip1j, hy_ip1j, tip1j);

					// check location and projections to determine force
					if ((tim1j < 0 || tim1j > 1) && (tip1j < 0 || tip1j > 1)){
						// IF vertex gj projects onto "bulk" of edge gi, no other edge
						ftmpi = evSoftAdhesionForce(gi, gj, hr_ij, hx_ij, hy_ij, tij, dfx, dfy);
						pi = gi;
						pxi = x[NDIM * gi] + tij * segX(gi, 0);
						pyi = x[NDIM * gi + 1] + tij * segX(gi, 1);
						delxi = hx_ij;
						delyi = hy_ij;	
						tiprint = tij;
						fxi = dfx;
						fyi = dfy;
					}				
					else if (tim1j > 0 && tim1j < 1) {
						if (hr_ij < hr_im1j){
							ftmpi = evSoftAdhesionForce(gi, gj, hr_ij, hx_ij, hy_ij, tij, dfx, dfy);		
							pi = gi;		
							pxi = x[NDIM * gi] + tij * segX(gi, 0);
							pyi = x[NDIM * gi + 1] + tij * segX(gi, 1);
							delxi = hx_ij;
							delyi = hy_ij;	
							tiprint = tij;
							fxi = dfx;
							fyi = dfy;
						}
						else{
							ftmpi = evSoftAdhesionForce(gim1, gj, hr_im1j, hx_im1j, hy_im1j, tim1j, dfx, dfy);		
							pi = gim1;
							pxi = x[NDIM * gim1] + tim1j * segX(gim1, 0);
							pyi = x[NDIM * gim1 + 1] + tim1j * segX(gim1, 1);
							delxi = hx_im1j;
							delyi = hy_im1j;
							tiprint = tim1j;
							fxi = dfx;
							fyi = dfy;
						}
					}
					else if (tip1j > 0 && tip1j < 1) {
						if (hr_ij < hr_ip1j){
							ftmpi = evSoftAdhesionForce(gi, gj, hr_ij, hx_ij, hy_ij, tij, dfx, dfy);	
							pi = gi;
							pxi = x[NDIM * gi] + tij * segX(gi, 0);
							pyi = x[NDIM * gi + 1] + tij * segX(gi, 1);
							delxi = hx_ij;
							delyi = hy_ij;		
							tiprint = tij;		
							fxi = dfx;
							fyi = dfy;
						}
						else {
							ftmpi = evSoftAdhesionForce(gip1, gj, hr_ip1j, hx_ip1j, hy_ip1j, tip1j, dfx, dfy);
							pi = gip1;
							pxi = x[NDIM * gip1] + tip1j * segX(gip1, 0);
							pyi = x[NDIM * gip1 + 1] + tip1j * segX(gip1, 1);
							delxi = hx_ip1j;
							delyi = hy_ip1j;
							tiprint = tip1j;
							fxi = dfx;
							fyi = dfy;
						}
					}
				}



				// -- force due to projection from vertex i onto edge j (Note, vertex-vertex overlap will already have been taken care of above)
				if (drijp1 > shellij && dr > shellij && hr_ji < shellij && tji > 0 && tji < 1){
					// compute projection onto adjacent edges
					hr_jm1i = edge2VertexDistance(gjm1, gi, hx_jm1i, hy_jm1i, tjm1i);
					hr_jp1i = edge2VertexDistance(gjp1, gi, hx_jp1i, hy_jp1i, tjp1i);

					// check location and projections to determine force
					if ((tjm1i < 0 || tjm1i > 1) && (tjp1i < 0 || tjp1i > 1)) {
						// IF vertex gi projects onto bulk of edge gi and no other edge
						ftmpj = evSoftAdhesionForce(gj, gi, hr_ji, hx_ji, hy_ji, tji, dfx, dfy);
						pj = gj;
						pxj = x[NDIM * gj] + tji * segX(gj, 0);
						pyj = x[NDIM * gj + 1] + tji * segX(gj, 1);
						delxj = hx_ji;
						delyj = hy_ji;
						tjprint = tji;
						fxj = dfx;
						fyj = dfy;
					}
					else if (tjm1i > 0 && tjm1i < 1){
						if (hr_ji < hr_jm1i){
							ftmpj = evSoftAdhesionForce(gj, gi, hr_ji, hx_ji, hy_ji, tji, dfx, dfy);				// IF vertex gi projects onto edges gj and gjm1, closer to gj (REGION III')
							pj = gj;
							pxj = x[NDIM * gj] + tji * segX(gj, 0);
							pyj = x[NDIM * gj + 1] + tji * segX(gj, 1);
							delxj = hx_ji;
							delyj = hy_ji;
							tjprint = tji;
							fxj = dfx;
							fyj = dfy;
						}
						else {
							ftmpj = evSoftAdhesionForce(gjm1, gi, hr_jm1i, hx_jm1i, hy_jm1i, tjm1i, dfx, dfy);		// IF vertex gi projects onto edges gj and gjm1, closer to gjm1 (REGION II')
							pj = gjm1;
							pxj = x[NDIM * gjm1] + tjm1i * segX(gjm1, 0);
							pyj = x[NDIM * gjm1 + 1] + tjm1i * segX(gjm1, 1);
							delxj = hx_jm1i;
							delyj = hy_jm1i;
							tjprint = tjm1i;
							fxj = dfx;
							fyj = dfy;
						}
					}
					else if (tjp1i > 0 && tjp1i < 1) {
						if (hr_ji < hr_jp1i){
							ftmpj = evSoftAdhesionForce(gj, gi, hr_ji, hx_ji, hy_ji, tji, dfx, dfy);				// IF vertex gi projects onto edges gj and gjm1, closer to gj (REGION II' centered on ip1)
							pj = gj;
							pxj = x[NDIM * gj] + tji * segX(gj, 0);
							pyj = x[NDIM * gj + 1] + tji * segX(gj, 1);
							delxj = hx_ji;
							delyj = hy_ji;
							tjprint = tji;
							fxj = dfx;
							fyj = dfy;
						}
						else {
							ftmpj = evSoftAdhesionForce(gjp1, gi, hr_jp1i, hx_jp1i, hy_jp1i, tjp1i, dfx, dfy);		// IF vertex gi projects onto edges gj and gjm1, closer to gjm1 centered on ip1)
							pj = gjp1;
							pxj = x[NDIM * gjp1] + tjp1i * segX(gjp1, 0);
							pyj = x[NDIM * gjp1 + 1] + tjp1i * segX(gjp1, 1);
							delxj = hx_jp1i;
							delyj = hy_jp1i;
							tjprint = tjp1i;
							fxj = dfx;
							fyj = dfy;
						}
					}
				}



				// print based on forces
				if (ftmpi > 0 || ftmpi < 0){
					// print out index
					outputObj << pi << "  ";
					outputObj << gj << "  ";
					outputObj << setprecision(8) << pxi << "  ";
					outputObj << setprecision(8) << pyi << "  ";
					outputObj << setprecision(8) << delxi << "  ";
					outputObj << setprecision(8) << delyi << "  ";
					outputObj << setprecision(8) << tiprint << "  ";
					outputObj << setprecision(8) << fxi << "  ";
					outputObj << setprecision(8) << fyi << "  ";

					// print out torque
					if (tiprint > 0){
						// force acting between edge i and vertex j
						tqi = x[NDIM * pi] * (1 - tiprint) * fyi - x[NDIM * pi + 1] * (1 - tiprint) * fxi;
						tqi += x[NDIM * ip1[pi]] * tiprint * fyi - x[NDIM * ip1[pi] + 1] * tiprint * fxi;
						tqj = (pxi + delxi) * (-fyi) - (pyi + delyi) * (-fxi);
					}
					else {
						// based on force acting between vertex i and vertex j
						tqi = pxi * fyi - pyi * fxi;
						tqj = (pxi + delxi) * (-fyi) - (pyi + delyi) * (-fxi);
					}
					outputObj << setprecision(8) << tqi << "  ";
					outputObj << setprecision(8) << tqj << "  ";

					// print out new line character
					outputObj << endl;
				}
				if (ftmpj > 0 || ftmpj < 0){
					// print out easy stuff
					outputObj << pj << "  ";
					outputObj << gi << "  ";
					outputObj << setprecision(8) << pxj << "  ";
					outputObj << setprecision(8) << pyj << "  ";
					outputObj << setprecision(8) << delxj << "  ";
					outputObj << setprecision(8) << delyj << "  ";
					outputObj << setprecision(8) << tjprint << "  ";
					outputObj << setprecision(8) << fxj << "  ";
					outputObj << setprecision(8) << fyj << "  ";

					// print out torque
					if (tjprint > 0){
						// force acting between edge i and vertex j
						tqj = pxj * fyj - pyj * fxj;
						tqi = (pxj + delxj) * (-fyj) - (pyj + delyj) * (-fxj);
					}
					else {
						// based on force acting between vertex i and vertex j
						tqj = x[NDIM * pj] * (1 - tjprint) * fyj - x[NDIM * pj + 1] * (1 - tjprint) * fxj;
						tqj += x[NDIM * ip1[pj]] * tjprint * fyj - x[NDIM * ip1[pj] + 1] * tjprint * fxj;
						tqi = (pxj + delxj) * (-fyj) - (pyj + delyj) * (-fxj);
					}
					outputObj << setprecision(8) << tqi << "  ";
					outputObj << setprecision(8) << tqj << "  ";


					// print out new line character
					outputObj << endl;
				}
			}
		}
	}

	// print to close frame
	fill(F.begin(), F.end(), 0.0);
	outputObj << "ENDFR" << endl;
}

