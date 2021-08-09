/*

	BASIC FUNCTION DEFINITIONS for DPM class

	Jack Treado, 04/10/21

*/

#include "dpm.h"
#include <functional>

// namespace
using namespace Eigen;
using namespace std;

/******************************

	C O N S T R U C T O R S  & 

		D E S T R U C T O R

*******************************/

// Constructor with only one input (dim)
dpm::dpm(int ndim) {
	// local variables
	int d, i;

	// print to console
	cout << "** Instantiating dpm object, ndim = " << ndim << endl;

	// main variables
	NDIM = ndim;
	NNN = 4;

	// set scalars to default values
	dt = 0.0;

	ka = 0.0;
	kl = 0.0;
	kb = 0.0;
	kc = 0.0;

	l1 = 0.0;
	l2 = 0.0;

	// default boundary variables
	L.resize(NDIM);
	pbc.resize(NDIM);
	for (d = 0; d < NDIM; d++) {
		L[d] = 1.0;
		pbc[d] = 1;
	}

	// macroscopic stress vector
	stress.resize(NDIM * (NDIM + 1) / 2);
	for (i = 0; i < NDIM * (NDIM + 1) / 2; i++)
		stress.at(i) = 0.0;

	// initialize nearest neighbor info
	NBX = -1;
}

// Main constructor
dpm::dpm(int n, int ndim, int seed) {
	// local variables
	int d, i;

	// print to console
	cout << "** Instantiating dpm object, NCELLS = " << n << ",  ndim = " << ndim << ", seed = " << seed << " ..." << endl;

	// main variables
	NCELLS = n;
	NDIM = ndim;
	NNN = 4;

	// set scalars to default values
	dt = 0.0;

	ka = 0.0;
	kl = 0.0;
	kb = 0.0;
	kc = 0.0;

	l1 = 0.0;
	l2 = 0.0;

	// default boundary variables
	L.resize(NDIM);
	pbc.resize(NDIM);
	for (d = 0; d < NDIM; d++) {
		L[d] = 1.0;
		pbc[d] = 1;
	}

	// preferred area for each cell
	a0.resize(NCELLS);

	// macroscopic stress vector
	stress.resize(NDIM * (NDIM + 1) / 2);
	for (i = 0; i < NDIM * (NDIM + 1) / 2; i++)
		stress.at(i) = 0.0;

	// contact network vector
	cij.resize(NCELLS * (NCELLS - 1) / 2);
	for (i = 0; i < NCELLS * (NCELLS - 1) / 2; i++)
		cij.at(i) = 0;

	// initialize nearest neighbor info
	NBX = -1;

	// seed random number generator
	srand48(seed);
}

// destructor
dpm::~dpm() {
	// clear all private vectors
	L.clear();
	pbc.clear();
	a0.clear();
	l0.clear();
	t0.clear();
	nv.clear();
	szList.clear();
	im1.clear();
	ip1.clear();
	r.clear();
	x.clear();
	v.clear();
	F.clear();
	stress.clear();
	sb.clear();
	lb.clear();
	for (int i = 0; i < NBX; i++)
		nn.at(i).clear();
	nn.clear();
	head.clear();
	last.clear();
	list.clear();

	if (posout.is_open())
		posout.close();
}

/******************************

	C E L L   S H A P E

	G E T T E R S

*******************************/

// get global vertex index gi given input cell index ci and vertex index vi
int dpm::gindex(int ci, int vi) {
	return szList[ci] + vi;
}

// get cell index ci and vertex index
void dpm::cindices(int &ci, int &vi, int gi) {
	for (int i = NCELLS - 1; i >= 0; i--) {
		if (gi >= szList[i]) {
			ci = i;
			vi = gi - szList[ci];
			break;
		}
	}
}

// get cell area
double dpm::area(int ci) {
	// local variables
	int vi, vip1, gi, gip1, nvtmp;
	double dx, dy, xi, yi, xip1, yip1, areaVal = 0.0;

	// initial position: vi = 0
	nvtmp = nv.at(ci);
	gi = gindex(ci, 0);
	xi = x[NDIM * gi];
	yi = x[NDIM * gi + 1];

	// loop over vertices of cell ci, get area by shoe-string method
	for (vi = 0; vi < nvtmp; vi++) {
		// next vertex
		gip1 = ip1[gi];
		gi++;

		// get positions (check minimum images)
		dx = x[NDIM * gip1] - xi;
		if (pbc[0])
			dx -= L[0] * round(dx / L[0]);
		xip1 = xi + dx;

		dy = x[NDIM * gip1 + 1] - yi;
		if (pbc[1])
			dy -= L[1] * round(dy / L[1]);
		yip1 = yi + dy;

		// increment area
		areaVal += xi * yip1 - xip1 * yi;

		// set next coordinates
		xi = xip1;
		yi = yip1;
	}
	areaVal *= 0.5;

	return abs(areaVal);
}

// get cell perimeter
double dpm::perimeter(int ci) {
	// local variables
	int vi, gi, gip1, nvtmp;
	double dx, dy, xi, yi, xip1, yip1, l, perimVal = 0.0;

	// initial position: vi = 0
	nvtmp = nv.at(ci);
	gi = gindex(ci, 0);
	xi = x[NDIM * gi];
	yi = x[NDIM * gi + 1];

	// loop over vertices of cell ci, get perimeter
	for (vi = 0; vi < nvtmp; vi++) {
		// next vertex
		gip1 = ip1[gi];
		gi++;

		// get positions (check minimum images)
		dx = x[NDIM * gip1] - xi;
		if (pbc[0])
			dx -= L[0] * round(dx / L[0]);
		xip1 = xi + dx;

		dy = x[NDIM * gip1 + 1] - yi;
		if (pbc[1])
			dy -= L[1] * round(dy / L[1]);
		yip1 = yi + dy;

		// compute segment length
		l = sqrt(dx * dx + dy * dy);

		// add to perimeter
		perimVal += l;

		// update coordinates
		xi = xip1;
		yi = yip1;
	}

	// return perimeter
	return perimVal;
}

// get cell center of mass position
void dpm::com2D(int ci, double &cx, double &cy) {
	// local variables
	int vi, gi, gip1, nvtmp;
	double dx, dy, xi, yi, xip1, yip1, l;

	// initial position: vi = 0
	nvtmp = nv.at(ci);
	gi = gindex(ci, 0);
	xi = x[NDIM * gi];
	yi = x[NDIM * gi + 1];

	// initialize center of mass coordinates
	cx = xi;
	cy = yi;

	// loop over vertices of cell ci, get perimeter
	for (vi = 0; vi < nvtmp - 1; vi++) {
		// next vertex
		gip1 = ip1.at(gi);
		gi++;

		// get positions (check minimum images)
		dx = x[NDIM * gip1] - xi;
		if (pbc[0])
			dx -= L[0] * round(dx / L[0]);
		xip1 = xi + dx;

		dy = x[NDIM * gip1 + 1] - yi;
		if (pbc[1])
			dy -= L[1] * round(dy / L[1]);
		yip1 = yi + dy;

		// add to center of mass
		cx += xip1;
		cy += yip1;

		// update coordinates
		xi = xip1;
		yi = yip1;
	}

	// take average to get com
	cx /= nvtmp;
	cy /= nvtmp;
}

// get configuration packing fraction
double dpm::vertexPackingFraction2D() {
	int ci;
	double val, boxV, areaSum = 0.0;

	// numerator
	for (ci = 0; ci < NCELLS; ci++)
		areaSum += area(ci) + 0.25 * PI * pow(2.0 * r.at(szList[ci]), 2.0) * (0.5 * nv.at(ci) - 1);

	// denominator
	boxV = L[0] * L[1];

	// return packing fraction
	val = areaSum / boxV;
	return val;
}

// get configuration "preferred" packing fraction
double dpm::vertexPreferredPackingFraction2D() {
	int ci;
	double val, boxV, areaSum = 0.0;

	// numerator
	for (ci = 0; ci < NCELLS; ci++)
		areaSum += a0[ci] + 0.25 * PI * pow(2.0 * r.at(szList[ci]), 2.0) * (0.5 * nv.at(ci) - 1);

	// denominator
	boxV = L[0] * L[1];

	// return packing fraction
	val = areaSum / boxV;
	return val;
}

// get vertex kinetic energy
double dpm::vertexKineticEnergy() {
	double K = 0;

	for (int i = 0; i < vertDOF; i++)
		K += v[i] * v[i];
	K *= 0.5;

	return K;
}

// get number of vertex-vertex contacts
int dpm::vvContacts() {
	int nvv = 0;

	for (int ci = 0; ci < NCELLS; ci++) {
		for (int cj = ci + 1; cj < NCELLS; cj++)
			nvv += cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
	}

	return nvv;
}

// get number of cell-cell contacts
int dpm::ccContacts() {
	int ncc = 0;

	for (int ci = 0; ci < NCELLS; ci++) {
		for (int cj = ci + 1; cj < NCELLS; cj++) {
			if (cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2] > 0)
				ncc++;
		}
	}

	return ncc;
}

/******************************

	I N I T I A L -

			I Z A T I O N

*******************************/

// initialize vertex indexing
void dpm::initializeVertexIndexing2D() {
	int gi, vi, vip1, vim1, ci;

	// check that vertDOF has been assigned
	if (NVTOT <= 0){
		cerr << "	** ERROR: in initializeVertexIndexing2D, NVTOT not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cerr << "	** ERROR: in initializeVertexIndexing2D, vertDOF not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cerr << "	** ERROR: in initializeVertexIndexing2D, nv vector not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
		exit(1);
	}

	// save list of adjacent vertices
	im1.resize(NVTOT);
	ip1.resize(NVTOT);
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

	// initialize vertex configuration vectors
	x.resize(vertDOF);
	v.resize(vertDOF);
	F.resize(vertDOF);
}

// initialize vertex shape parameters based on nv
void dpm::initializeVertexShapeParameters(double calA0, int nref) {
	// local variables
	int gi, ci, vi, nvtmp;
	double lenscale, calA0tmp, calAntmp;

	// check that vertDOF has been assigned
	if (NVTOT <= 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, NVTOT not assigned. Ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, vertDOF not assigned. Ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, nv vector not assigned. Ending here." << endl;
    	exit(1);
	}

	// resize shape paramters
	l0.resize(NVTOT);
	t0.resize(NVTOT);
	r.resize(NVTOT);

	// loop over cells, give each cell the same shape parameters based on calA0 + nv
	for (ci = 0; ci < NCELLS; ci++){
		lenscale = (double) nv.at(ci)/nref;
		initializeVertexShapeParameters(ci,calA0,lenscale);
	}
}


// initialize vertex shape parameters for SPECIFIC CELL ci
void dpm::initializeVertexShapeParameters(int ci, double calA0, double lenscale) {
	// local variables
	int gi, vi, nvtmp;
	double calA0tmp, calAntmp;

	// check that vertDOF has been assigned
	if (NVTOT <= 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, NVTOT not assigned. Ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, vertDOF not assigned. Ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cerr << "	** ERROR: in initializeVertexShapeParameters, nv vector not assigned. Ending here." << endl;
    	exit(1);
	}

	// check that ci is within bounds
	if (ci >= NCELLS){
		cerr << "	** ERROR: in initializeVertexShapeParameters, ci = " << ci << ", but NCELLS = " << NCELLS << ". Will cause out of bounds error, so ending here. " << endl;
		exit(1);
	}

	// check that l0, t0 and r have been resized correctly
	if (l0.size() != NVTOT){
		cerr << "	** ERROR: in initializeVertexShapeParameters, l0 vector assigned incorrectly. l0.size() = " << l0.size() << ", while NVTOT = " << NVTOT << ". Ending here." << endl;
    	exit(1);
	}
	if (t0.size() != NVTOT){
		cerr << "	** ERROR: in initializeVertexShapeParameters, t0 vector assigned incorrectly. t0.size() = " << t0.size() << ", while NVTOT = " << NVTOT << ". Ending here." << endl;
    	exit(1);
	}
	if (r.size() != NVTOT){
		cerr << "	** ERROR: in initializeVertexShapeParameters, r vector assigned incorrectly. r.size() = " << r.size() << ", while NVTOT = " << NVTOT << ". Ending here." << endl;
    	exit(1);
	}

	// since everything has been pre-allocated at this point, determine shape parameters based on input calA0
	nvtmp = nv.at(ci);

	// a0 based on input lengthscale
	a0.at(ci) = lenscale * lenscale;

	// shape parameter
	calAntmp = nvtmp * tan(PI / nvtmp) / PI;
	calA0tmp = calA0 * calAntmp;

	// l0 and vertex radii
	gi = szList.at(ci);
	for (vi = 0; vi < nvtmp; vi++) {
		l0.at(gi + vi) = 2.0 * sqrt(PI * calA0tmp * a0.at(ci)) / nvtmp;
		t0.at(gi + vi) = 0.0;
		r.at(gi + vi) = 0.5 * l0.at(gi + vi);
	}
}

// initialize bidisperse cell system, single calA0
void dpm::bidisperse2D(double calA0, int nsmall, double smallfrac, double sizefrac) {
	// local variables
	double calA0tmp, calAntmp, rtmp, areaSum;
	int vim1, vip1, gi, ci, vi, nlarge, smallN, largeN, NVSMALL;

	// print to console
	cout << "** initializing bidisperse DPM particles in 2D ..." << endl;

	// number of vertices on large particles
	nlarge = round(sizefrac * nsmall);

	// total number of vertices
	smallN = round(smallfrac * NCELLS);
	largeN = NCELLS - smallN;
	NVSMALL = nsmall * smallN;
	NVTOT = NVSMALL + nlarge * largeN;
	vertDOF = NDIM * NVTOT;

	// szList and nv (keep track of global vertex indices)
	nv.resize(NCELLS);
	szList.resize(NCELLS);

	nv.at(0) = nsmall;
	for (ci = 1; ci < NCELLS; ci++){
		if (ci < smallN) {
			nv.at(ci) = nsmall;
			szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
		}
		else {
			nv.at(ci) = nlarge;
			szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);
		}
	}

	// initialize vertex shape parameters
	initializeVertexShapeParameters(calA0, nsmall);

	// initialize vertex indexing
	initializeVertexIndexing2D();
}

// initialize gaussian polydisperse cell system, single calA0
void dpm::gaussian2D(double dispersion, double calA0, int n1) {
	// local variables
	double calA0tmp, calAntmp, rtmp, areaSum, r1, r2, grv;
	int vim1, vip1, gi, ci, vi, nvtmp;

	// print to console
	cout << "** initializing gaussian DPM particles in 2D with size dispersion " << dispersion << " ..." << endl;

	// szList and nv (keep track of global vertex indices)
	nv.resize(NCELLS);
	szList.resize(NCELLS);

	nv.at(0) = n1;
	NVTOT = n1;
	for (ci = 1; ci < NCELLS; ci++) {
		// use Box-Muller to generate polydisperse sample
		r1 = drand48();
		r2 = drand48();
		grv = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
		nvtmp = floor(dispersion * n1 * grv + n1);
		if (nvtmp < nvmin)
			nvtmp = nvmin;

		// store size of cell ci
		nv.at(ci) = nvtmp;
		szList.at(ci) = szList.at(ci - 1) + nv.at(ci - 1);

		// add to total NV count
		NVTOT += nvtmp;
	}
	vertDOF = NDIM * NVTOT;

	// initialize vertex shape parameters (give all cells same calA0 / calAn)
	initializeVertexShapeParameters(calA0, n1);

	// initialize vertex indexing
	initializeVertexIndexing2D();
}

// set sinusoidal preferred angle
void dpm::sinusoidalPreferredAngle(double thA, double thK) {
	int ci, vi, gi;
	double thR;

	// print to console
	cout << "** setting initial th0 values to sinusoids, thA = " << thA << ", thK = " << thK << " ..." << endl;

	// loop over cells
	gi = 0;
	for (ci = 0; ci < NCELLS; ci++) {
		thR = (2.0 * PI) / nv.at(ci);
		for (vi = 0; vi < nv.at(ci); vi++) {
			t0.at(gi) = thA * thR * sin(thR * thK * vi);
			gi++;
		}
	}
}

// initialize CoM positions of cells using SP FIRE
void dpm::initializePositions2D(double phi0, double Ftol) {
	// local variables
	int i, d, ci, cj, vi, vj, gi, cellDOF = NDIM * NCELLS;
	double areaSum, xtra = 1.1;

	// local disk vectors
	vector<double> drad(NCELLS, 0.0);
	vector<double> dpos(cellDOF, 0.0);
	vector<double> dv(cellDOF, 0.0);
	vector<double> dF(cellDOF, 0.0);

	// print to console
	cout << "** initializing particle positions using 2D SP model and FIRE relaxation ..." << endl;

	// initialize box size based on packing fraction
	areaSum = 0.0;
	for (ci = 0; ci < NCELLS; ci++)
		areaSum += a0.at(ci) + 0.25 * PI * pow(l0.at(ci), 2.0) * (0.5 * nv.at(ci) - 1);

	// set box size
	for (d = 0; d < NDIM; d++)
		L.at(d) = pow(areaSum / phi0, 1.0 / NDIM);

	// initialize cell centers randomly
	for (ci = 0; ci < cellDOF; ci += 2)
		dpos.at(ci) = L[ci % 2] * drand48();
	for (ci = cellDOF - 1; ci > 0; ci -= 2)
		dpos.at(ci) = L[ci % 2] * drand48();

	// set radii of SP disks
	for (ci = 0; ci < NCELLS; ci++)
		drad.at(ci) = xtra * sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin(2.0 * PI / nv.at(ci))));

	// FIRE VARIABLES
	double P = 0;
	double fnorm = 0;
	double vnorm = 0;
	double alpha = alpha0;

	double dt0 = 1e-2;
	double dtmax = 10 * dt0;
	double dtmin = 1e-8 * dt0;

	int npPos = 0;
	int npNeg = 0;

	int fireit = 0;
	double fcheck = 10 * Ftol;

	// interaction variables
	double rij, sij, dtmp, ftmp, vftmp;
	double dr[NDIM];

	// initial step size
	dt = dt0;

	// loop until force relaxes
	while ((fcheck > Ftol) && fireit < itmax) {
		// FIRE step 1. Compute P
		P = 0.0;
		for (i = 0; i < cellDOF; i++)
			P += dv[i] * dF[i];

		// FIRE step 2. adjust simulation based on net motion of degrees of freedom
		if (P > 0) {
			// increase positive counter
			npPos++;

			// reset negative counter
			npNeg = 0;

			// alter simulation if enough positive steps have been taken
			if (npPos > NMIN) {
				// change time step
				if (dt * finc < dtmax)
					dt *= finc;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else {
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
			for (i = 0; i < cellDOF; i++)
			{
				// take half step backwards
				dpos[i] -= 0.5 * dt * dv[i];

				// reset velocities
				dv[i] = 0.0;
			}

			// decrease time step if past initial delay
			if (fireit > NDELAY)
			{
				// decrease time step
				if (dt * fdec > dtmin)
					dt *= fdec;

				// reset alpha
				alpha = alpha0;
			}
		}

		// FIRE step 3. First VV update
		for (i = 0; i < cellDOF; i++)
			dv[i] += 0.5 * dt * dF[i];

		// FIRE step 4. adjust velocity magnitude
		fnorm = 0.0;
		vnorm = 0.0;
		for (i = 0; i < cellDOF; i++) {
			fnorm += dF[i] * dF[i];
			vnorm += dv[i] * dv[i];
		}
		fnorm = sqrt(fnorm);
		vnorm = sqrt(vnorm);
		if (fnorm > 0) {
			for (i = 0; i < cellDOF; i++)
				dv[i] = (1 - alpha) * dv[i] + alpha * (vnorm / fnorm) * dF[i];
		}

		// FIRE step 4. Second VV update
		for (i = 0; i < cellDOF; i++) {
			dpos[i] += dt * dv[i];
			dF[i] = 0.0;
		}

		// FIRE step 5. Update forces
		for (ci = 0; ci < NCELLS; ci++) {
			for (cj = ci + 1; cj < NCELLS; cj++) {

				// contact distance
				sij = drad[ci] + drad[cj];

				// true distance
				rij = 0.0;
				for (d = 0; d < NDIM; d++) {
					// get distance element
					dtmp = dpos[NDIM * cj + d] - dpos[NDIM * ci + d];
					if (pbc[d])
						dtmp -= L[d] * round(dtmp / L[d]);

					// add to true distance
					rij += dtmp * dtmp;

					// save in distance array
					dr[d] = dtmp;
				}
				rij = sqrt(rij);

				// check distances
				if (rij < sij) {
					// force magnitude
					ftmp = kc * (1.0 - (rij / sij)) / sij;

					// add to vectorial force
					for (d = 0; d < NDIM; d++)
					{
						vftmp = ftmp * (dr[d] / rij);
						dF[NDIM * ci + d] -= vftmp;
						dF[NDIM * cj + d] += vftmp;
					}
				}
			}
		}

		// FIRE step 5. Final VV update
		for (i = 0; i < cellDOF; i++)
			dv[i] += 0.5 * dt * dF[i];

		// update forces to check
		fcheck = 0.0;
		for (i = 0; i < cellDOF; i++)
			fcheck += dF[i] * dF[i];
		fcheck = sqrt(fcheck / NCELLS);

		// print to console
		if (fireit % NSKIP == 0) {
			cout << endl
				 << endl;
			cout << "===========================================" << endl;
			cout << "		I N I T I A L  S P 			" << endl;
			cout << " 	F I R E 						" << endl;
			cout << "		M I N I M I Z A T I O N 	" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit = " << fireit << endl;
			cout << "	** fcheck = " << fcheck << endl;
			cout << "	** fnorm = " << fnorm << endl;
			cout << "	** vnorm = " << vnorm << endl;
			cout << "	** dt = " << dt << endl;
			cout << "	** P = " << P << endl;
			cout << "	** Pdir = " << P / (fnorm * vnorm) << endl;
			cout << "	** alpha = " << alpha << endl;
		}

		// update iterate
		fireit++;
	}
	// check if FIRE converged
	if (fireit == itmax) {
		cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
		exit(1);
	}
	else {
		cout << endl
			 << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===========================================" << endl;
		cout << " 	F I R E 						" << endl;
		cout << "		M I N I M I Z A T I O N 	" << endl;
		cout << "	C O N V E R G E D! 				" << endl
			 << endl;

		cout << "	(for initial disk minimization) " << endl;
		cout << "===========================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	** fireit = " << fireit << endl;
		cout << "	** fcheck = " << fcheck << endl;
		cout << "	** vnorm = " << vnorm << endl;
		cout << "	** dt = " << dt << endl;
		cout << "	** P = " << P << endl;
		cout << "	** alpha = " << alpha << endl;
	}

	// initialize vertex positions based on cell centers
	for (ci = 0; ci < NCELLS; ci++) {
		for (vi = 0; vi < nv.at(ci); vi++) {
			// get global vertex index
			gi = gindex(ci, vi);

			// length from center to vertex
			dtmp = sqrt((2.0 * a0.at(ci)) / (nv.at(ci) * sin((2.0 * PI) / nv.at(ci))));

			// set positions
			x.at(NDIM * gi) = dtmp * cos((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci) + 1e-2 * l0[gi] * drand48();
			x.at(NDIM * gi + 1) = dtmp * sin((2.0 * PI * vi) / nv.at(ci)) + dpos.at(NDIM * ci + 1) + 1e-2 * l0[gi] * drand48();
		}
	}
}

// initialize neighbor linked list
void dpm::initializeNeighborLinkedList2D(double boxLengthScale) {
	// local variables
	double llscale;
	int i, d, nntmp, scx;

	// print to console
	cout << "** initializing neighbor linked list, boxLengthScale = " << boxLengthScale;

	// get largest radius as llscale
	llscale = 2.0 * (*max_element(r.begin(),r.end()));

	// initialize box length vectors
	NBX = 1;
	sb.resize(NDIM);
	lb.resize(NDIM);
	for (d = 0; d < NDIM; d++) {
		// determine number of cells along given dimension by rmax
		sb[d] = round(L[d] / (boxLengthScale * llscale));

		// just in case, if < 3, change to 3 so box neighbor checking will work
		if (sb[d] < 3)
			sb[d] = 3;

		// determine box length by number of cells
		lb[d] = L[d] / sb[d];

		// count total number of cells
		NBX *= sb[d];
	}

	// initialize list of box nearest neighbors
	scx = sb[0];
	nn.resize(NBX);

	// loop over cells, save forward neighbors for each box
	for (i = 0; i < NBX; i++) {
		// reshape entry
		nn[i].resize(NNN);

		// right neighbor (i+1)
		nn[i][0] = (i + 1) % NBX; 

		// top neighbors (i,j+1), (i+1,j+1)
		if (pbc[1]){
			// (i,j+1) w/ pbc
			nn[i][1] = (i + scx) % NBX;

			// (i+1,j+1) w/ pbc
			nn[i][2] = (nn[i][1] + 1) % NBX;
		}
		else {
			// if on top row, both = -1
			if (i >= NBX - scx){
				nn[i][1] = -1;
				nn[i][2] = -1;
			}
			// if not on top row, still add
			else{
				nn[i][1] = i + scx; 
				nn[i][2] = nn[i][1] + 1;
			}
		}

		// bottom neighbor w/ pbc (j-1)
		nntmp = (i + NBX - scx) % NBX;	

		// bottom-right neighbor (i+1, j-1)
		if (pbc[1])
			nn[i][3] = nntmp + 1;
		else{
			// if on bottom row, skip
			if (i < scx)
				nn[i][3] = -1;
			// otherwise, set
			else
				nn[i][3] = nntmp + 1;
		}

		// right-hand bc (periodic)
		if ((i + 1) % scx == 0) {
			if (pbc[0]) {
				nn[i][0] = i - scx + 1;
				if (pbc[1]){
					nn[i][2] = nn[i][1] - scx + 1;
					nn[i][3] = nntmp - scx + 1;
				}
			}
			else {
				nn[i][0] = -1;
				nn[i][2] = -1;
				nn[i][3] = -1;
			}
		}
	}

	// linked-list variables
	head.resize(NBX);
	last.resize(NBX);
	list.resize(NVTOT + 1);

	// print box info to console
	cout << ";  initially NBX = " << NBX << " ..." << endl;
}

/******************************

	E D I T I N G   &

			U P D A T I N G

*******************************/

// sort vertices into neighbor linked list
void dpm::sortNeighborLinkedList2D() {
	// local variables
	int d, gi, boxid, sbtmp;
	double xtmp;

	// reset linked list info
	fill(list.begin(), list.end(), 0);
	fill(head.begin(), head.end(), 0);
	fill(last.begin(), last.end(), 0);

	// sort vertices into linked list
	for (gi = 0; gi < NVTOT; gi++) {
		// 1. get cell id of current particle position
		boxid = 0;
		sbtmp = 1;
		for (d = 0; d < NDIM; d++) {
			// current location
			xtmp = x[NDIM * gi + d];

			// check out-of-bounds
			if (xtmp < 0){
				if (pbc[d])
					xtmp += L[d];
				else
					xtmp = 0.00001;
			}
			else if (xtmp > L[d]){
				if (pbc[d])
					xtmp -= L[d];
				else
					xtmp = 0.99999*L[d];
			}

			// add d index to 1d list
			boxid += floor(xtmp / lb[d]) * sbtmp;

			// increment dimensional factor
			sbtmp *= sb[d];
		}

		// 2. add to head list or link within list
		// NOTE: particle ids are labelled starting from 1, setting to 0 means end of linked list
		if (head[boxid] == 0) {
			head[boxid] = gi + 1;
			last[boxid] = gi + 1;
		}
		else {
			list[last[boxid]] = gi + 1;
			last[boxid] = gi + 1;
		}
	}
}

// change size of particles
void dpm::scaleParticleSizes2D(double scaleFactor) {
	// local variables
	int gi, ci, vi, xind, yind;
	double xi, yi, cx, cy, dx, dy;

	// loop over cells, scale
	for (ci = 0; ci < NCELLS; ci++) {
		// scale preferred area
		a0[ci] *= scaleFactor * scaleFactor;

		// first global index for ci
		gi = szList.at(ci);

		// compute cell center of mass
		xi = x[NDIM * gi];
		yi = x[NDIM * gi + 1];
		cx = xi;
		cy = yi;
		for (vi = 1; vi < nv.at(ci); vi++) {
			dx = x.at(NDIM * (gi + vi)) - xi;
			if (pbc[0])
				dx -= L[0] * round(dx / L[0]);

			dy = x.at(NDIM * (gi + vi) + 1) - yi;
			if (pbc[1])
				dy -= L[1] * round(dy / L[1]);

			xi += dx;
			yi += dy;

			cx += xi;
			cy += yi;
		}
		cx /= nv.at(ci);
		cy /= nv.at(ci);

		for (vi = 0; vi < nv.at(ci); vi++) {
			// x and y inds
			xind = NDIM * (gi + vi);
			yind = xind + 1;

			// closest relative position
			dx = x[xind] - cx;
			if (pbc[0])
				dx -= L[0] * round(dx / L[0]);

			dy = x[yind] - cy;
			if (pbc[1])
				dy -= L[1] * round(dy / L[1]);

			// update vertex positions
			x[xind] += (scaleFactor - 1.0) * dx;
			x[yind] += (scaleFactor - 1.0) * dy;

			// scale vertex radii
			r[gi + vi] *= scaleFactor;
			l0[gi + vi] *= scaleFactor;
		}
	}
}

// remove rattlers from contact network, return number of rattlers
int dpm::removeRattlers() {
	// local variables
	int ci, cj, ctmp, rvv, rcc, nr, nm = 1;

	// loop over rows, eliminate contacts to rattlers until nm = 0
	while (nm > 0) {
		// reset number of rattlers
		nr = 0;

		// number of "marginal" rattlers to be removed
		nm = 0;
		for (ci = 0; ci < NCELLS; ci++) {
			// get number of contacts on cell ci
			rvv = 0;
			rcc = 0;
			for (cj = 0; cj < NCELLS; cj++) {
				if (ci != cj) {
					if (ci > cj)
						ctmp = cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];
					else
						ctmp = cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
				}
				else
					ctmp = 0;

				rvv += ctmp;
				if (ctmp > 0)
					rcc++;
			}

			// check to see if particle should be removed from network
			if (rcc <= NDIM && rvv <= 3) {
				// increment # of rattlers
				nr++;

				// if in contact, remove contacts
				if (rvv > 0) {
					nm++;

					for (cj = 0; cj < NCELLS; cj++) {
						// delete contact between ci and cj
						if (ci != cj) {
							if (ci > cj)
								cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2] = 0;
							else
								cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2] = 0;
						}
					}
				}
			}
		}
	}

	// return total number of rattlers
	return nr;
}

// draw random velocities based on input temperature
void dpm::drawVelocities2D(double T) {
	// local variables
	int gi;
	double r1, r2, grv1, grv2, gnorm, tscale = sqrt(T), vcomx = 0.0, vcomy = 0.0;

	// loop over velocities, draw from maxwell-boltzmann distribution
	for (gi = 0; gi < NVTOT; gi++) {
		// draw random numbers using Box-Muller
		r1 = drand48();
		r2 = drand48();
		grv1 = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
		grv2 = sqrt(-2.0 * log(r1)) * sin(2.0 * PI * r2);

		// gaussian velocities norm
		gnorm = sqrt(grv1*grv1 + grv2*grv2);

		// assign to velocities
		v[NDIM * gi] = tscale * (grv1/gnorm);
		v[NDIM * gi + 1] = tscale * (grv2/gnorm);

		// add to center of mass
		vcomx += v[NDIM * gi];
		vcomy += v[NDIM * gi + 1];
	}
	vcomx = vcomx / NVTOT;
	vcomy = vcomy / NVTOT;

	// subtract off center of mass drift
	for (gi = 0; gi < NVTOT; gi++) {
		v[NDIM * gi] -= vcomx;
		v[NDIM * gi + 1] -= vcomy;
	}
}

/******************************

	D P M  F O R C E 

			U P D A T E S

*******************************/

void dpm::resetForcesAndEnergy() {
	fill(F.begin(), F.end(), 0.0);
	fill(stress.begin(), stress.end(), 0.0);
	U = 0.0;
}

void dpm::shapeForces2D() {
	// local variables
	int ci, gi, vi, nvtmp;
	double fa, fli, flim1, fb, cx, cy, xi, yi;
	double rho0, l0im1, l0i, a0tmp, atmp;
	double dx, dy, da, dli, dlim1, dtim1, dti, dtip1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
	double nim1x, nim1y, nix, niy, sinim1, sini, sinip1, cosim1, cosi, cosip1;
	double ddtim1, ddti;

	// loop over vertices, add to force
	rho0 = sqrt(a0.at(0));
	ci = 0;
	for (gi = 0; gi < NVTOT; gi++) {

		// -- Area force (and get cell index ci)
		if (ci < NCELLS) {
			if (gi == szList[ci]) {
				// shape information
				nvtmp = nv[ci];
				a0tmp = a0[ci];

				// preferred segment length of last segment
				l0im1 = l0[im1[gi]];

				// compute area deviation
				atmp = area(ci);
				da = (atmp / a0tmp) - 1.0;

				// update potential energy
				U += 0.5 * ka * (da * da);

				// shape force parameters
				fa = ka * da * (rho0 / a0tmp);
				fb = kb * rho0;

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
				rix = x[NDIM * gi] - cx;
				riy = x[NDIM * gi + 1] - cy;

				// get prior adjacent vertices
				rim2x = x[NDIM * im1[im1[gi]]] - cx;
				rim2y = x[NDIM * im1[im1[gi]] + 1] - cy;
				if (pbc[0])
					rim2x -= L[0] * round(rim2x / L[0]);
				if (pbc[1])
					rim2y -= L[1] * round(rim2y / L[1]);

				rim1x = x[NDIM * im1[gi]] - cx;
				rim1y = x[NDIM * im1[gi] + 1] - cy;
				if (pbc[0])
					rim1x -= L[0] * round(rim1x / L[0]);
				if (pbc[1])
					rim1y -= L[1] * round(rim1y / L[1]);

				// get prior segment vectors
				lim2x = rim1x - rim2x;
				lim2y = rim1y - rim2y;

				lim1x = rix - rim1x;
				lim1y = riy - rim1y;

				// increment cell index
				ci++;
			}
		}

		// preferred segment length
		l0i = l0[gi];

		// get next adjacent vertices
		rip1x = x[NDIM * ip1[gi]] - cx;
		rip1y = x[NDIM * ip1[gi] + 1] - cy;
		if (pbc[0])
			rip1x -= L[0] * round(rip1x / L[0]);
		if (pbc[1])
			rip1y -= L[1] * round(rip1y / L[1]);

		// -- Area force
		F[NDIM * gi] += 0.5 * fa * (rim1y - rip1y);
		F[NDIM * gi + 1] += 0.5 * fa * (rip1x - rim1x);

		// -- Perimeter force
		lix = rip1x - rix;
		liy = rip1y - riy;

		// segment lengths
		lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);
		li = sqrt(lix * lix + liy * liy);

		// segment deviations
		dlim1 = (lim1 / l0im1) - 1.0;
		dli = (li / l0i) - 1.0;

		// segment forces
		flim1 = kl * (rho0 / l0im1);
		fli = kl * (rho0 / l0i);

		// add to forces
		F[NDIM * gi] += (fli * dli * lix / li) - (flim1 * dlim1 * lim1x / lim1);
		F[NDIM * gi + 1] += (fli * dli * liy / li) - (flim1 * dlim1 * lim1y / lim1);

		// update potential energy
		U += 0.5 * kl * (dli * dli);

		// -- Bending force
		if (kb > 0) {
			// get ip2 for third angle
			rip2x = x[NDIM * ip1[ip1[gi]]] - cx;
			rip2y = x[NDIM * ip1[ip1[gi]] + 1] - cy;
			if (pbc[0])
				rip2x -= L[0] * round(rip2x / L[0]);
			if (pbc[1])
				rip2y -= L[1] * round(rip2y / L[1]);

			// get last segment length
			lip1x = rip2x - rip1x;
			lip1y = rip2y - rip1y;

			// get angles
			sinim1 = lim1x * lim2y - lim1y * lim2x;
			cosim1 = lim1x * lim2x + lim1y * lim2y;

			sini = lix * lim1y - liy * lim1x;
			cosi = lix * lim1x + liy * lim1y;

			sinip1 = lip1x * liy - lip1y * lix;
			cosip1 = lip1x * lix + lip1y * liy;

			// get normal vectors
			nim1x = lim1y;
			nim1y = -lim1x;

			nix = liy;
			niy = -lix;

			// get change in angles
			dtim1 = atan2(sinim1, cosim1) - t0[im1[gi]];
			dti = atan2(sini, cosi) - t0[gi];
			dtip1 = atan2(sinip1, cosip1) - t0[ip1[gi]];

			// get delta delta theta's
			ddtim1 = (dti - dtim1) / (lim1 * lim1);
			ddti = (dti - dtip1) / (li * li);

			// add to force
			F[NDIM * gi] += fb * (ddtim1 * nim1x + ddti * nix);
			F[NDIM * gi + 1] += fb * (ddtim1 * nim1y + ddti * niy);

			// update potential energy
			U += 0.5 * kb * (dti * dti);
		}

		// update old coordinates
		rim2x = rim1x;
		rim1x = rix;
		rix = rip1x;

		rim2y = rim1y;
		rim1y = riy;
		riy = rip1y;

		// update old segment vectors
		lim2x = lim1x;
		lim2y = lim1y;

		lim1x = lix;
		lim1y = liy;

		// update old preferred segment length
		l0im1 = l0i;
	}
}

void dpm::vertexRepulsiveForces2D() {
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, dx, dy, rho0;
	double ftmp, fx, fy;

	// sort particles
	sortNeighborLinkedList2D();

	// get fundamental length
	rho0 = sqrt(a0.at(0));

	// reset contact network
	fill(cij.begin(), cij.end(), 0);

	// loop over boxes in neighbor linked list
	for (bi = 0; bi < NBX; bi++) {
		// get start of list of vertices
		pi = head[bi];

		// loop over linked list
		while (pi > 0) {
			// real particle index
			gi = pi - 1;

			// next particle in list
			pj = list[pi];

			// loop down neighbors of pi in same cell
			while (pj > 0) {
				// real index of pj
				gj = pj - 1;

				if (gj == ip1[gi] || gj == im1[gi]) {
					pj = list[pj];
					continue;
				}

				// contact distance
				sij = r[gi] + r[gj];

				// particle distance
				dx = x[NDIM * gj] - x[NDIM * gi];
				if (pbc[0])
					dx -= L[0] * round(dx / L[0]);
				if (dx < sij) {
					dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
					if (pbc[1])
						dy -= L[1] * round(dy / L[1]);
					if (dy < sij) {
						rij = sqrt(dx * dx + dy * dy);
						if (rij < sij) {
							// force scale
							ftmp = kc * (1 - (rij / sij)) * (rho0 / sij);
							fx = ftmp * (dx / rij);
							fy = ftmp * (dy / rij);

							// add to forces
							F[NDIM * gi] -= fx;
							F[NDIM * gi + 1] -= fy;

							F[NDIM * gj] += fx;
							F[NDIM * gj + 1] += fy;

							// increae potential energy
							U += 0.5 * kc * pow((1 - (rij / sij)), 2.0);

							// add to virial stress
							stress[0] += dx * fx;
							stress[1] += dy * fy;
							stress[2] += 0.5 * (dx * fy + dy * fx);

							// add to contacts
							cindices(ci, vi, gi);
							cindices(cj, vj, gj);

							if (ci > cj)
								cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
							else if (ci < cj)
								cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
						}
					}
				}

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
					// contact distance
					sij = r[gi] + r[gj];

					// particle distance
					dx = x[NDIM * gj] - x[NDIM * gi];
					if (pbc[0])
						dx -= L[0] * round(dx / L[0]);
					if (dx < sij) {
						dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
						if (pbc[1])
							dy -= L[1] * round(dy / L[1]);
						if (dy < sij) {
							rij = sqrt(dx * dx + dy * dy);
							if (rij < sij) {
								// force scale
								ftmp = kc * (1 - (rij / sij)) * (rho0 / sij);
								fx = ftmp * (dx / rij);
								fy = ftmp * (dy / rij);

								// add to forces
								F[NDIM * gi] -= fx;
								F[NDIM * gi + 1] -= fy;

								F[NDIM * gj] += fx;
								F[NDIM * gj + 1] += fy;

								// increae potential energy
								U += 0.5 * kc * pow((1 - (rij / sij)), 2.0);

								// add to virial stress
								stress[0] += dx * fx;
								stress[1] += dy * fy;
								stress[2] += 0.5 * (dx * fy + dy * fx);

								// add to contacts
								cindices(ci, vi, gi);
								cindices(cj, vj, gj);

								if (ci > cj)
									cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2]++;
								else if (ci < cj)
									cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;
							}
						}
					}

					// update pj
					pj = list[pj];
				}
			}

			// update pi index to be next
			pi = list[pi];
		}
	}

	// normalize stress by box area, make dimensionless
	stress[0] *= (rho0 / (L[0] * L[1]));
	stress[1] *= (rho0 / (L[0] * L[1]));
	stress[2] *= (rho0 / (L[0] * L[1]));
}

void dpm::vertexAttractiveForces2D() {
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, dx, dy, rho0;
	double ftmp, fx, fy;

	// attraction shell parameters
	double shellij, cutij, xij, kint = (kc*l1)/(l2 - l1);

	// sort particles
	sortNeighborLinkedList2D();

	// get fundamental length
	rho0 = sqrt(a0[0]);

	// reset contact network
	fill(cij.begin(), cij.end(), 0);

	// loop over boxes in neighbor linked list
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

				if (gj == ip1[gi] || gj == im1[gi]) {
					pj = list[pj];
					continue;
				}

				// cell index of gj
				cindices(cj, vj, gj);

				// contact distance
				sij = r[gi] + r[gj];

				// attraction distances
				shellij = (1.0 + l2)*sij;
				cutij = (1.0 + l1)*sij;

				// particle distance
				dx = x[NDIM * gj] - x[NDIM * gi];
				if (pbc[0])
					dx -= L[0] * round(dx / L[0]);
				if (dx < shellij) {
					dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
					if (pbc[1])
						dy -= L[1] * round(dy / L[1]);
					if (dy < shellij) {
						rij = sqrt(dx * dx + dy * dy);
						if (rij < shellij) {
							// scaled distance
							xij = rij/sij;

							// pick force based on vertex-vertex distance
							if (rij > cutij){
								// force scale
								ftmp = kint*(xij - 1.0 - l2)/sij;

								// increase potential energy
								U += -0.5*kint*pow(1.0 + l2 - xij,2.0);
							}
							else{
								// force scale
								ftmp = kc*(1 - xij)/sij;

								// increase potential energy
								U += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
							}

							// force elements
							fx 					= ftmp*(dx/rij);
							fy 					= ftmp*(dy/rij);

							// add to forces
							F[NDIM*gi] 			-= fx;
							F[NDIM*gi + 1] 		-= fy;

							F[NDIM*gj] 			+= fx;
							F[NDIM*gj + 1] 		+= fy;

							// add to virial stress
							stress[0] 			+= dx*fx;
							stress[1] 			+= dy*fy;
							stress[2] 			+= 0.5*(dx*fy + dy*fx);

							// add to contacts
							if (ci > cj)
								cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
							else if (ci < cj)
								cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++;
						}
					}
				}

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

					// cell index of gj
					cindices(cj, vj, gj);

					// contact distance
					sij = r[gi] + r[gj];

					// attraction distances
					shellij = (1.0 + l2)*sij;
					cutij = (1.0 + l1)*sij;

					// particle distance
					dx = x[NDIM * gj] - x[NDIM * gi];
					if (pbc[0])
						dx -= L[0] * round(dx / L[0]);
					if (dx < shellij) {
						dy = x[NDIM * gj + 1] - x[NDIM * gi + 1];
						if (pbc[1])
							dy -= L[1] * round(dy / L[1]);
						if (dy < shellij) {
							rij = sqrt(dx * dx + dy * dy);
							if (rij < shellij) {
								// scaled distance
								xij = rij/sij;

								// pick force based on vertex-vertex distance
								if (rij > cutij){
									// force scale
									ftmp = kint*(xij - 1.0 - l2)/sij;

									// increase potential energy
									U += -0.5*kint*pow(1.0 + l2 - xij,2.0);
								}
								else{
									// force scale
									ftmp = kc*(1 - xij)/sij;

									// increase potential energy
									U += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
								}

								// force elements
								fx 					= ftmp*(dx/rij);
								fy 					= ftmp*(dy/rij);

								// add to forces
								F[NDIM*gi] 			-= fx;
								F[NDIM*gi + 1] 		-= fy;

								F[NDIM*gj] 			+= fx;
								F[NDIM*gj + 1] 		+= fy;

								// add to virial stress
								stress[0] 			+= dx*fx;
								stress[1] 			+= dy*fy;
								stress[2] 			+= 0.5*(dx*fy + dy*fx);

								if (ci > cj)
									cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
								else if (ci < cj)
									cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++; 
							}
						}
					}

					// update pj
					pj = list[pj];
				}
			}

			// update pi index to be next
			pi = list[pi];
		}
	}

	// normalize stress by box area, make dimensionless
	stress[0] *= (rho0 / (L[0] * L[1]));
	stress[1] *= (rho0 / (L[0] * L[1]));
	stress[2] *= (rho0 / (L[0] * L[1]));
}





void dpm::repulsiveForceUpdate() {
	resetForcesAndEnergy();
	shapeForces2D();
	vertexRepulsiveForces2D();
}


void dpm::attractiveForceUpdate(){
	resetForcesAndEnergy();
	shapeForces2D();
	vertexAttractiveForces2D();
}



/******************************

	D P M  

		I N T E G R A T O R S

*******************************/

void dpm::setdt(double dt0) {
	// local variables
	int i;
	double ta, tl, tb, tmin, rho0;

	// typical length
	rho0 = sqrt(a0.at(0));

	// set typical time scales
	ta = rho0 / sqrt(ka);
	tl = (rho0 * l0.at(0)) / sqrt(ka * kl);
	tb = (rho0 * l0.at(0)) / sqrt(ka * kb);

	// set main time scale as min
	tmin = 1e8;
	if (ta < tmin)
		tmin = ta;
	if (tl < tmin)
		tmin = tl;
	if (tb < tmin)
		tmin = tb;

	// set dt
	dt = dt0 * tmin;
}

void dpm::vertexFIRE2D(dpmMemFn forceCall, double Ftol, double dt0) {
	// local variables
	int i;
	double rho0;

	// check to see if cell linked-list has been initialized
	if (NBX == -1) {
		cerr << "	** ERROR: In dpm::fire, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
		exit(1);
	}

	// FIRE variables
	double P, fnorm, fcheck, vnorm, alpha, dtmax, dtmin;
	int npPos, npNeg, fireit;

	// set dt based on geometric parameters
	setdt(dt0);

	// Initialize FIRE variables
	P = 0;
	fnorm = 0;
	vnorm = 0;
	alpha = alpha0;

	dtmax = 10.0 * dt;
	dtmin = 1e-2 * dt;

	npPos = 0;
	npNeg = 0;

	fireit = 0;
	fcheck = 10 * Ftol;

	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);

	// length scale
	rho0 = sqrt(a0.at(0));

	// relax forces using FIRE
	while (fcheck > Ftol && fireit < itmax) {
		// compute P
		P = 0.0;
		for (i = 0; i < vertDOF; i++)
			P += v[i] * F[i];

		// print to console
		if (fireit % NSKIP == 0) {
			cout << endl
				 << endl;
			cout << "===========================================" << endl;
			cout << " 	F I R E 						" << endl;
			cout << "		M I N I M I Z A T I O N 	" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit 	= " << fireit << endl;
			cout << "	** fcheck 	= " << fcheck << endl;
			cout << "	** U 		= " << U << endl;
			cout << "	** dt 		= " << dt << endl;
			cout << "	** P 		= " << P << endl;
			cout << "	** alpha 	= " << alpha << endl;
			cout << "	** npPos 	= " << npPos << endl;
			cout << "	** npNeg 	= " << npNeg << endl;
			cout << "	** sxx  	= " << stress[0] << endl;
			cout << "	** syy 		= " << stress[1] << endl;
			cout << "	** sxy 		= " << stress[2] << endl;
		}

		// Adjust simulation based on net motion of degrees of freedom
		if (P > 0) {
			// increase positive counter
			npPos++;

			// reset negative counter
			npNeg = 0;

			// alter simulation if enough positive steps have been taken
			if (npPos > NDELAY) {
				// change time step
				if (dt * finc < dtmax)
					dt *= finc;

				// decrease alpha
				alpha *= falpha;
			}
		}
		else {
			// reset positive counter
			npPos = 0;

			// increase negative counter
			npNeg++;

			// check if simulation is stuck
			if (npNeg > NNEGMAX) {
				cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
				exit(1);
			}

			// take half step backwards, reset velocities
			for (i = 0; i < vertDOF; i++) {
				// take half step backwards
				x[i] -= 0.5 * dt * v[i];

				// reset vertex velocities
				v[i] = 0.0;
			}

			// decrease time step if past initial delay
			if (fireit > NDELAY) {
				// decrease time step
				if (dt * fdec > dtmin)
					dt *= fdec;

				// reset alpha
				alpha = alpha0;
			}
		}

		// VV VELOCITY UPDATE #1
		for (i = 0; i < vertDOF; i++)
			v[i] += 0.5 * dt * F[i];

		// compute fnorm, vnorm and P
		fnorm = 0.0;
		vnorm = 0.0;
		for (i = 0; i < vertDOF; i++) {
			fnorm += F[i] * F[i];
			vnorm += v[i] * v[i];
		}
		fnorm = sqrt(fnorm);
		vnorm = sqrt(vnorm);

		// update velocities (s.d. vs inertial dynamics) only if forces are acting
		if (fnorm > 0) {
			for (i = 0; i < vertDOF; i++)
				v[i] = (1 - alpha) * v[i] + alpha * (F[i] / fnorm) * vnorm;
		}

		// VV POSITION UPDATE
		for (i = 0; i < vertDOF; i++) {
			// update position
			x[i] += dt * v[i];

			// recenter in box
			if (x[i] > L[i % NDIM] && pbc[i % NDIM])
				x[i] -= L[i % NDIM];
			else if (x[i] < 0 && pbc[i % NDIM])
				x[i] += L[i % NDIM];
		}

		// update forces (function passed as argument)
		CALL_MEMBER_FN(*this, forceCall)();

		// VV VELOCITY UPDATE #2
		for (i = 0; i < vertDOF; i++)
			v[i] += 0.5 * F[i] * dt;

		// update fcheck based on fnorm (= force per degree of freedom)
		fcheck = 0.0;
		for (i = 0; i < vertDOF; i++)
			fcheck += F[i] * F[i];
		fcheck = sqrt(fcheck / vertDOF);

		// update iterator
		fireit++;
	}
	// check if FIRE converged
	if (fireit == itmax) {
		cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
		exit(1);
	}
	else {
		cout << endl;
		cout << "===========================================" << endl;
		cout << " 	F I R E 						" << endl;
		cout << "		M I N I M I Z A T I O N 	" << endl;
		cout << "	C O N V E R G E D! 				" << endl;
		cout << "===========================================" << endl;
		cout << endl;
		cout << "	** fireit 	= " << fireit << endl;
		cout << "	** fcheck 	= " << fcheck << endl;
		cout << "	** U 		= " << U << endl;

		cout << "	** fnorm	= " << fnorm << endl;
		cout << "	** vnorm 	= " << vnorm << endl;
		cout << "	** dt 		= " << dt << endl;
		cout << "	** P 		= " << P << endl;
		cout << "	** alpha 	= " << alpha << endl;
		cout << "	** sxx  	= " << stress[0] << endl;
		cout << "	** syy 		= " << stress[1] << endl;
		cout << "	** sxy 		= " << stress[2] << endl;
		cout << endl << endl;
	}
}

void dpm::vertexNVE2D(ofstream &enout, dpmMemFn forceCall, double T, double dt0, int NT, int NPRINTSKIP){
	// local variables
	int t, i;
	double K, simclock;

	// set time step magnitude
	setdt(dt0);

	// initialize time keeper
	simclock = 0.0;

	// initialize velocities
	drawVelocities2D(T);

	// loop over time, print energy
	for (t=0; t<NT; t++){
		// VV VELOCITY UPDATE #1
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*dt*F[i];

		// VV POSITION UPDATE
		for (i=0; i<vertDOF; i++){
			// update position
			x[i] += dt*v[i];

			// recenter in box
			if (x[i] > L[i % NDIM] && pbc[i % NDIM])
				x[i] -= L[i % NDIM];
			else if (x[i] < 0 && pbc[i % NDIM])
				x[i] += L[i % NDIM];
		}

		// FORCE UPDATE
		CALL_MEMBER_FN(*this, forceCall)();

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*F[i]*dt;

		// update sim clock
		simclock += dt;

		// print to console and file
		if (t % NPRINTSKIP == 0){
			// compute kinetic energy
			K = vertexKineticEnergy();

			// print to console
			cout << endl << endl;
			cout << "===============================" << endl;
			cout << "	D P M  						" << endl;
			cout << " 			 					" << endl;
			cout << "		N V E 					" << endl;
			cout << "===============================" << endl;
			cout << endl;
			cout << "	** t / NT	= " << t << " / " << NT << endl;
			cout << "	** U 		= " << setprecision(12) << U << endl;
			cout << "	** K 		= " << setprecision(12) << K << endl;
			cout << "	** E 		= " << setprecision(12) << U + K << endl;

			// print to energy file
			cout << "** printing energy" << endl;
			enout << setw(w) << left << t;
			enout << setw(wnum) << left << simclock;
			enout << setw(wnum) << setprecision(12) << U;
			enout << setw(wnum) << setprecision(12) << K;
			enout << endl;

			// print to configuration only if position file is open
			if (posout.is_open())
				printConfiguration2D();
		}
	}
}

// Langevin dynamics for printing
void dpm::vertexLangevinNVT2D(ofstream &enout, dpmMemFn forceCall, double T0, double gam, double dt0, int NT, int NPRINTSKIP){
	// local variables
	int t, i;
	double T, K, simclock, dmp1, dmp2;
	double r1, r2, grv1, grv2;

	// set time step magnitude
	setdt(dt0);

	// set damping coefficients
	dmp1 = exp(-gam*dt);
	dmp2 = sqrt(1.0 - exp(-2.0*gam*dt))*sqrt(T0);

	// initialize time keeper
	simclock = 0.0;

	// initialize velocities
	drawVelocities2D(T0);

	// loop over time, print energy
	for (t=0; t<NT; t++){
		// Langevin velocity update #1
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*dt*F[i];

		// Langevin position update #1
		for (i=0; i<vertDOF; i++){
			// update position
			x[i] += 0.5*dt*v[i];

			// recenter in box
			if (x[i] > L[i % NDIM] && pbc[i % NDIM])
				x[i] -= L[i % NDIM];
			else if (x[i] < 0 && pbc[i % NDIM])
				x[i] += L[i % NDIM];
		}

		// Langevin random velocity update
		// use Box-Muller to generate gaussian random variables
		for (i=0; i<NVTOT; i++){
			r1 = drand48();
			r2 = drand48();
			grv1 = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
			grv2 = sqrt(-2.0 * log(r1)) * sin(2.0 * PI * r2);

			// add to both dimensions
			v[NDIM*i] = dmp1*v[NDIM*i] + dmp2*grv1;
			v[NDIM*i + 1] = dmp1*v[NDIM*i + 1] + dmp2*grv2;
		}

		// Langevin position update #2 (based on random kick)
		for (i=0; i<vertDOF; i++){
			// update position
			x[i] += 0.5*dt*v[i];

			// recenter in box
			if (x[i] > L[i % NDIM] && pbc[i % NDIM])
				x[i] -= L[i % NDIM];
			else if (x[i] < 0 && pbc[i % NDIM])
				x[i] += L[i % NDIM];
		}

		// FORCE UPDATE
		CALL_MEMBER_FN(*this, forceCall)();

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*F[i]*dt;

		// update sim clock
		simclock += dt;

		// print to console and file
		if (t % NPRINTSKIP == 0){
			// compute kinetic energy
			K = vertexKineticEnergy();

			// kinetic temperature
			T = K/NVTOT;

			// print to console
			cout << endl << endl;
			cout << "===============================" << endl;
			cout << "	D P M  						" << endl;
			cout << " 			 					" << endl;
			cout << "	L A N G E V I N 			" << endl;
			cout << "===============================" << endl;
			cout << endl;
			cout << "	** t / NT	= " << t << " / " << NT << endl;
			cout << "	** U 		= " << setprecision(12) << U << endl;
			cout << "	** K 		= " << setprecision(12) << K << endl;
			cout << "	** E 		= " << setprecision(12) << U + K << endl;
			cout << "	** T  		= " << setprecision(12) << T << endl;

			// print to energy file
			cout << "** printing energy" << endl;
			enout << setw(w) << left << t;
			enout << setw(wnum) << left << simclock;
			enout << setw(wnum) << setprecision(12) << U;
			enout << setw(wnum) << setprecision(12) << K;
			enout << setw(wnum) << setprecision(12) << T;
			enout << endl;

			// print to configuration only if position file is open
			if (posout.is_open())
				printConfiguration2D();
		}
	}
}

// Langevin dynamics with no printing
void dpm::vertexLangevinNVT2D(dpmMemFn forceCall, double T0, double gam, double dt0, int NT, int NPRINTSKIP){
	// local variables
	int t, i;
	double T, K, simclock, dmp1, dmp2;
	double r1, r2, grv1, grv2;

	// set time step magnitude
	setdt(dt0);

	// set damping coefficients
	dmp1 = exp(-gam*dt);
	dmp2 = sqrt(1.0 - exp(-2.0*gam*dt))*sqrt(T0);

	// initialize time keeper
	simclock = 0.0;

	// initialize velocities
	drawVelocities2D(T0);

	// loop over time, print energy
	for (t=0; t<NT; t++){
		// Langevin velocity update #1
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*dt*F[i];

		// Langevin position update #1
		for (i=0; i<vertDOF; i++){
			// update position
			x[i] += 0.5*dt*v[i];

			// recenter in box
			if (x[i] > L[i % NDIM] && pbc[i % NDIM])
				x[i] -= L[i % NDIM];
			else if (x[i] < 0 && pbc[i % NDIM])
				x[i] += L[i % NDIM];
		}

		// Langevin random velocity update
		// use Box-Muller to generate gaussian random variables
		for (i=0; i<NVTOT; i++){
			r1 = drand48();
			r2 = drand48();
			grv1 = sqrt(-2.0 * log(r1)) * cos(2.0 * PI * r2);
			grv2 = sqrt(-2.0 * log(r1)) * sin(2.0 * PI * r2);

			// add to both dimensions
			v[NDIM*i] = dmp1*v[NDIM*i] + dmp2*grv1;
			v[NDIM*i + 1] = dmp1*v[NDIM*i + 1] + dmp2*grv2;
		}

		// Langevin position update #2 (based on random kick)
		for (i=0; i<vertDOF; i++){
			// update position
			x[i] += 0.5*dt*v[i];

			// recenter in box
			if (x[i] > L[i % NDIM] && pbc[i % NDIM])
				x[i] -= L[i % NDIM];
			else if (x[i] < 0 && pbc[i % NDIM])
				x[i] += L[i % NDIM];
		}

		// FORCE UPDATE
		CALL_MEMBER_FN(*this, forceCall)();

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*F[i]*dt;

		// update sim clock
		simclock += dt;

		// print to console and file
		if (t % NPRINTSKIP == 0){
			// compute kinetic energy
			K = vertexKineticEnergy();

			// kinetic temperature
			T = K/NVTOT;

			// print to console
			cout << endl << endl;
			cout << "===============================" << endl;
			cout << "	D P M  						" << endl;
			cout << " 			 					" << endl;
			cout << "	L A N G E V I N 			" << endl;
			cout << "===============================" << endl;
			cout << endl;
			cout << "	** t / NT	= " << t << " / " << NT << endl;
			cout << "	** U 		= " << setprecision(12) << U << endl;
			cout << "	** K 		= " << setprecision(12) << K << endl;
			cout << "	** E 		= " << setprecision(12) << U + K << endl;
			cout << "	** T  		= " << setprecision(12) << T << endl;
		}
	}
}


/******************************

	D P M  

		P R O T O C O L S

*******************************/

void dpm::vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0) {
	// local variables
	int it = 0, itmax = 1e4;
	double phi0 = vertexPreferredPackingFraction2D();
	double scaleFactor=1.0, P, Sxy;

	// loop while phi0 < phi0Target
	while (phi0 < phi0Target && it < itmax) {
		// scale particle sizes
		scaleParticleSizes2D(scaleFactor);

		// update phi0
		phi0 = vertexPreferredPackingFraction2D();

		// relax configuration (pass member function force update)
		vertexFIRE2D(forceCall, Ftol, dt0);

		// get scale factor
		scaleFactor = sqrt((phi0 + dphi0) / phi0);

		// get updated pressure
		P = 0.5 * (stress[0] + stress[1]);
		Sxy = stress[2];

		// print to console
		cout << endl
			 << endl;
		cout << "===============================" << endl;
		cout << "								" << endl;
		cout << " 	C O M P R E S S I O N 		" << endl;
		cout << "								" << endl;
		cout << "	P R O T O C O L 	  		" << endl;
		cout << "								" << endl;
		cout << "===============================" << endl;
		cout << endl;
		cout << "	** it 			= " << it << endl;
		cout << "	** phi0 curr	= " << phi0 << endl;
		if (phi0 + dphi0 < phi0Target)
			cout << "	** phi0 next 	= " << phi0 + dphi0 << endl;
		cout << "	** P 			= " << P << endl;
		cout << "	** Sxy 			= " << Sxy << endl;
		cout << "	** U 			= " << U << endl;
		printConfiguration2D();
		cout << endl
			 << endl;

		// update iterate
		it++;
	}
}

void dpm::vertexJamming2D(dpmMemFn forceCall, double Ftol, double Ptol, double dt0, double dphi0, bool plotCompression) {
	// local variables
	int k = 0, nr;
	bool jammed = 0, overcompressed = 0, undercompressed = 0, stalled = 0;
	double pcheck, phi0, rH, r0, rL, rho0, scaleFactor = 1.0;

	// initialize binary root search parameters
	r0 = sqrt(a0.at(0));
	rH = -1;
	rL = -1;

	// initialize preferred packing fraction
	phi0 = vertexPreferredPackingFraction2D();

	// save initial state
	vector<double> xsave(vertDOF, 0.0);
	vector<double> vsave(vertDOF, 0.0);
	vector<double> Fsave(vertDOF, 0.0);
	vector<double> rsave(vertDOF, 0.0);
	vector<double> l0save(vertDOF, 0.0);
	vector<double> t0save(vertDOF, 0.0);
	vector<double> a0save(vertDOF, 0.0);

	xsave = x;
	rsave = r;
	l0save = l0;
	t0save = t0;
	a0save = a0;

	// loop until jamming is found
	while (!jammed && k < itmax) {
		// set length scale by 1st particle preferred area
		rho0 = sqrt(a0.at(0));

		// relax configuration (pass member function force update)
		vertexFIRE2D(forceCall, Ftol, dt0);

		// update pressure
		pcheck = 0.5 * (stress[0] + stress[1]);

		// remove rattlers
		nr = removeRattlers();

		// boolean checks for jamming
		undercompressed = ((pcheck < 2.0 * Ptol && rH < 0) || (pcheck < Ptol && rH > 0));
		overcompressed = (pcheck > 2.0 * Ptol);
		jammed = (pcheck < 2.0 * Ptol && pcheck > Ptol && rH > 0 && rL > 0);
		stalled = (rH > 0 && rL > 0 && abs(rH - rL) < 1e-8);

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << " 	Q U A S I S T A T I C  						" << endl;
		cout << " 	  	I S O T R O P I C 						" << endl;
		cout << "			C O M P R E S S I O N 				" << endl << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi0 			= " << phi0 << endl;
		cout << "	* phi 			= " << vertexPackingFraction2D() << endl;
		cout << "	* scaleFactor 	= " << scaleFactor << endl;
		cout << "	* r0 			= " << r0 << endl;
		cout << "	* rH 			= " << rH << endl;
		cout << "	* rL 			= " << rL << endl;
		cout << "	* pcheck 		= " << pcheck << endl;
		cout << "	* U 		 	= " << U << endl;
		cout << "	* Nvv  			= " << vvContacts() << endl;
		cout << "	* Ncc 			= " << ccContacts() << endl;
		cout << "	* # of rattlers = " << nr << endl << endl;
		cout << "	* undercompressed = " << undercompressed << endl;
		cout << "	* overcompressed = " << overcompressed << endl;
		cout << "	* jammed = " << jammed << endl;
		cout << "	* stalled = " << stalled << endl << endl;
		if (plotCompression)
			printConfiguration2D();
		cout << endl << endl;

		// update particle scaleFactor based on target check
		if (rH < 0) {
			// if still undercompressed, then grow until overcompressed found
			if (undercompressed) {
				r0 = rho0;
				scaleFactor = sqrt((phi0 + dphi0) / phi0);
			}
			// if first overcompressed, decompress by dphi/2 until unjamming
			else if (overcompressed) {
				// current = upper bound length scale r
				rH = rho0;

				// save first overcompressed state
				r0 = rH;
				xsave = x;
				rsave = r;
				vsave = v;
				Fsave = F;
				l0save = l0;
				t0save = t0;
				a0save = a0;

				// shrink particle sizes
				scaleFactor = sqrt((phi0 - 0.5 * dphi0) / phi0);

				// print to console
				cout << "	-- -- overcompressed for the first time, scaleFactor = " << scaleFactor << endl;
			}
		}
		else {
			if (rL < 0) {
				// if first undercompressed, save last overcompressed state, begin root search
				if (undercompressed) {
					// current = new lower bound length scale r
					rL = rho0;

					// load state
					x = xsave;
					r = rsave;
					v = vsave;
					F = Fsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute new scale factor by root search
					scaleFactor = 0.5 * (rH + rL) / r0;

					// print to console
					cout << "	-- -- undercompressed for the first time, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- BEGINNING ROOT SEARCH IN ENTHALPY MIN PROTOCOL..." << endl;
					cout << "	-- -- rH = " << rH << endl;
					cout << " 	-- -- rL = " << rL << endl;
					cout << "	-- -- |rH - rL| = " << abs(rH - rL) << endl;
				}
				// if still overcompressed, decrement again
				else if (overcompressed) {
					// current = upper bound length scale r
					rH = rho0;

					// save overcompressed state
					r0 = rH;
					xsave = x;
					rsave = r;
					vsave = v;
					Fsave = F;
					l0save = l0;
					t0save = t0;
					a0save = a0;

					// keep shrinking at same rate until unjamming
					scaleFactor = sqrt((phi0 - 0.5 * dphi0) / phi0);

					// print to console
					cout << "	-- -- overcompressed, still no unjamming, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- rH = " << rH << endl;
					cout << " 	-- -- rL = " << rL << endl;
					cout << "	-- -- |rH - rL| = " << abs(rH - rL) << endl;
				}
			}
			else {
				if (stalled) {
					cout << "Simulation STALLED ... resetting to initial overcompression ... " << endl;
					cout << "\t * rH = " << rH << endl;
					cout << "\t * rL = " << rL << endl;
					cout << "\t * |rH - rL| = " << abs(rH - rL) << endl;

					// reset 
					rH = -1;
					rL = -1;

					// load state
					x = xsave;
					r = rsave;
					v = vsave;
					F = Fsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute phi0
					phi0 = vertexPreferredPackingFraction2D();

					// reset by decompression
					scaleFactor = sqrt((phi0 - 0.25 * dphi0) / phi0);
				}
				// if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
				else if (undercompressed) {
					// current = new lower bound length scale r
					rL = rho0;

					// load state
					x = xsave;
					r = rsave;
					v = vsave;
					F = Fsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute new scale factor by root search
					scaleFactor = 0.5 * (rH + rL) / r0;

					// print to console
					cout << "	-- -- undercompressed, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- rH = " << rH << endl;
					cout << " 	-- -- rL = " << rL << endl;
					cout << "	-- -- |rH - rL| = " << abs(rH - rL) << endl;
				}
				else if (overcompressed) {
					// current = upper bound length scale r
					rH = rho0;

					// load state
					x = xsave;
					r = rsave;
					v = vsave;
					F = Fsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute new scale factor
					scaleFactor = 0.5 * (rH + rL) / r0;

					// print to console
					cout << "	-- -- overcompressed, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- rH = " << rH << endl;
					cout << " 	-- -- rL = " << rL << endl;
					cout << "	-- -- |rH - rL| = " << abs(rH - rL) << endl;
				}
				else if (jammed) {
					cout << "	** At k = " << k << ", target pressure found!" << endl;
					cout << " WRITING ENTHALPY-MINIMIZED CONFIG TO FILE" << endl;
					cout << " ENDING COMPRESSION SIMULATION" << endl;
					scaleFactor = 1.0;
					break;
				}
			}
		}

		// scale particle sizes
		scaleParticleSizes2D(scaleFactor);

		// update packing fraction
		phi0 = vertexPreferredPackingFraction2D();

		// update iterate
		k++;
	}
}


//  add cooling protocol (determine delTemp given cooling to Tmin over # = coolsteps)
void dpm::vertexAnneal2Jam2D(dpmMemFn forceCall, double Ftol, double Ptol, double dt0, double dphi0, double T0, double trun, bool plotCompression){
	// local variables
	int k = 0, nr, NT = 1, NPRINTSKIP = 1, NTCOOL = 1, coolsteps = 10, coolt;
	bool jammed = 0, overcompressed = 0, undercompressed = 0, stalled = 0;
	double pcheck, phi0, rH, r0, rL, rho0, scaleFactor = 1.0;
	double gam = 1.0;
	double Tcool = T0, Tmin = 1e-8, delTemp;

	// initialize binary root search parameters
	r0 = sqrt(a0.at(0));
	rH = -1;
	rL = -1;

	// initialize preferred packing fraction
	phi0 = vertexPreferredPackingFraction2D();

	// save initial state
	vector<double> xsave(vertDOF, 0.0);
	vector<double> vsave(vertDOF, 0.0);
	vector<double> Fsave(vertDOF, 0.0);
	vector<double> rsave(vertDOF, 0.0);
	vector<double> l0save(vertDOF, 0.0);
	vector<double> t0save(vertDOF, 0.0);
	vector<double> a0save(vertDOF, 0.0);

	xsave = x;
	rsave = r;
	l0save = l0;
	t0save = t0;
	a0save = a0;

	// loop until jamming is found
	while (!jammed && k < itmax) {
		// set length scale by 1st particle preferred area
		rho0 = sqrt(a0.at(0));

		// run Langevin dynamics only if not root searching
		if (rH < 0){
			// get dt
			setdt(dt0);

			// run NVT for trun time at temperature T
			NT = (int) floor(trun/dt);
			NPRINTSKIP = (int) floor((0.05*trun)/dt);
			vertexLangevinNVT2D(forceCall, T0, gam, dt0, NT, NPRINTSKIP);			

			// cool down to T = 1e-8
			NTCOOL = NT/coolsteps;
			delTemp = pow(10.0,(log10(Tmin) - log10(T0))/coolsteps);
			Tcool = T0;
			for (coolt=0; coolt<coolsteps; coolt++){
				// equilibrate at lower temperature
				Tcool *= delTemp;
				cout << "** Equilibrating now at T = " << Tcool << endl;
				vertexLangevinNVT2D(forceCall, Tcool, gam, dt0, NTCOOL, NPRINTSKIP);
			}
		}
		
		// relax configuration (pass member function force update)
		vertexFIRE2D(forceCall, Ftol, dt0);

		// update pressure
		pcheck = 0.5 * (stress[0] + stress[1]);

		// remove rattlers
		nr = removeRattlers();

		// boolean checks for jamming
		undercompressed = ((pcheck < 2.0 * Ptol && rH < 0) || (pcheck < Ptol && rH > 0));
		overcompressed = (pcheck > 2.0 * Ptol);
		jammed = (pcheck < 2.0 * Ptol && pcheck > Ptol && rH > 0 && rL > 0);
		stalled = (rH > 0 && rL > 0 && abs(rH - rL) < 1e-8);

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << " 	Q U A S I S T A T I C  						" << endl;
		cout << " 	  	I S O T R O P I C 						" << endl;
		cout << "			C O M P R E S S I O N 				" << endl << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi0 			= " << phi0 << endl;
		cout << "	* phi 			= " << vertexPackingFraction2D() << endl;
		cout << "	* scaleFactor 	= " << scaleFactor << endl;
		cout << "	* r0 			= " << r0 << endl;
		cout << "	* rH 			= " << rH << endl;
		cout << "	* rL 			= " << rL << endl;
		cout << "	* pcheck 		= " << pcheck << endl;
		cout << "	* U 		 	= " << U << endl;
		cout << "	* Nvv  			= " << vvContacts() << endl;
		cout << "	* Ncc 			= " << ccContacts() << endl;
		cout << "	* # of rattlers = " << nr << endl << endl;
		cout << "	* undercompressed = " << undercompressed << endl;
		cout << "	* overcompressed = " << overcompressed << endl;
		cout << "	* jammed = " << jammed << endl;
		cout << "	* stalled = " << stalled << endl << endl;
		if (plotCompression)
			printConfiguration2D();
		cout << endl << endl;

		// update particle scaleFactor based on target check
		if (rH < 0) {
			// if still undercompressed, then grow until overcompressed found
			if (undercompressed) {
				r0 = rho0;
				scaleFactor = sqrt((phi0 + dphi0) / phi0);
			}
			// if first overcompressed, decompress by dphi/2 until unjamming
			else if (overcompressed) {
				// current = upper bound length scale r
				rH = rho0;

				// save first overcompressed state
				r0 = rH;
				xsave = x;
				rsave = r;
				vsave = v;
				Fsave = F;
				l0save = l0;
				t0save = t0;
				a0save = a0;

				// shrink particle sizes
				scaleFactor = sqrt((phi0 - 0.5 * dphi0) / phi0);

				// print to console
				cout << "	-- -- overcompressed for the first time, scaleFactor = " << scaleFactor << endl;
			}
		}
		else {
			if (rL < 0) {
				// if first undercompressed, save last overcompressed state, begin root search
				if (undercompressed) {
					// current = new lower bound length scale r
					rL = rho0;

					// load state
					x = xsave;
					r = rsave;
					v = vsave;
					F = Fsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute new scale factor by root search
					scaleFactor = 0.5 * (rH + rL) / r0;

					// print to console
					cout << "	-- -- undercompressed for the first time, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- BEGINNING ROOT SEARCH IN ENTHALPY MIN PROTOCOL..." << endl;
					cout << "	-- -- rH = " << rH << endl;
					cout << " 	-- -- rL = " << rL << endl;
					cout << "	-- -- |rH - rL| = " << abs(rH - rL) << endl;
				}
				// if still overcompressed, decrement again
				else if (overcompressed) {
					// current = upper bound length scale r
					rH = rho0;

					// save overcompressed state
					r0 = rH;
					xsave = x;
					rsave = r;
					vsave = v;
					Fsave = F;
					l0save = l0;
					t0save = t0;
					a0save = a0;

					// keep shrinking at same rate until unjamming
					scaleFactor = sqrt((phi0 - 0.5 * dphi0) / phi0);

					// print to console
					cout << "	-- -- overcompressed, still no unjamming, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- rH = " << rH << endl;
					cout << " 	-- -- rL = " << rL << endl;
					cout << "	-- -- |rH - rL| = " << abs(rH - rL) << endl;
				}
			}
			else {
				if (stalled) {
					cout << "Simulation STALLED ... resetting to initial overcompression ... " << endl;
					cout << "\t * rH = " << rH << endl;
					cout << "\t * rL = " << rL << endl;
					cout << "\t * |rH - rL| = " << abs(rH - rL) << endl;

					// reset 
					rH = -1;
					rL = -1;

					// load state
					x = xsave;
					r = rsave;
					v = vsave;
					F = Fsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute phi0
					phi0 = vertexPreferredPackingFraction2D();

					// reset by decompression
					scaleFactor = sqrt((phi0 - 0.25 * dphi0) / phi0);
				}
				// if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
				else if (undercompressed) {
					// current = new lower bound length scale r
					rL = rho0;

					// load state
					x = xsave;
					r = rsave;
					v = vsave;
					F = Fsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute new scale factor by root search
					scaleFactor = 0.5 * (rH + rL) / r0;

					// print to console
					cout << "	-- -- undercompressed, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- rH = " << rH << endl;
					cout << " 	-- -- rL = " << rL << endl;
					cout << "	-- -- |rH - rL| = " << abs(rH - rL) << endl;
				}
				else if (overcompressed) {
					// current = upper bound length scale r
					rH = rho0;

					// load state
					x = xsave;
					r = rsave;
					v = vsave;
					F = Fsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute new scale factor
					scaleFactor = 0.5 * (rH + rL) / r0;

					// print to console
					cout << "	-- -- overcompressed, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- rH = " << rH << endl;
					cout << " 	-- -- rL = " << rL << endl;
					cout << "	-- -- |rH - rL| = " << abs(rH - rL) << endl;
				}
				else if (jammed) {
					cout << "	** At k = " << k << ", target pressure found!" << endl;
					cout << " WRITING ENTHALPY-MINIMIZED CONFIG TO FILE" << endl;
					cout << " ENDING COMPRESSION SIMULATION" << endl;
					scaleFactor = 1.0;
					break;
				}
			}
		}

		// scale particle sizes
		scaleParticleSizes2D(scaleFactor);

		// update packing fraction
		phi0 = vertexPreferredPackingFraction2D();

		// update iterate
		k++;
	}
}

/******************************

	D P M  

		H E S S I A N

*******************************/

// wrapper function for total hessian
// note: dynamical matrix M = H - S
void dpm::dpmHessian2D(Eigen::MatrixXd &H, Eigen::MatrixXd &S) {
	// local variables
	int k, l;

	// print something to the console
	cout << "** Computing Hessian for configuration in dpmHessian2D ..." << endl;

	// initialize all possible matrices
	Eigen::MatrixXd Ha(vertDOF, vertDOF);  // stiffness matrix for area term
	Eigen::MatrixXd Sa(vertDOF, vertDOF);  // stress matrix for area term
	Eigen::MatrixXd Hl(vertDOF, vertDOF);  // stiffness matrix for perimeter term
	Eigen::MatrixXd Sl(vertDOF, vertDOF);  // stress matrix for perimeter term
	Eigen::MatrixXd Hb(vertDOF, vertDOF);  // stiffness matrix for bending energy
	Eigen::MatrixXd Sb(vertDOF, vertDOF);  // stress matrix for bending term
	Eigen::MatrixXd Hvv(vertDOF, vertDOF); // stiffness matrix for interaction terms
	Eigen::MatrixXd Svv(vertDOF, vertDOF); // stress matrix for interaction terms

	// initialize all matrices to be 0 initially
	for (k = 0; k < vertDOF; k++) {
		for (l = 0; l < vertDOF; l++) {
			Ha(k, l) = 0.0;
			Sa(k, l) = 0.0;
			Hl(k, l) = 0.0;
			Sl(k, l) = 0.0;
			Hb(k, l) = 0.0;
			Sb(k, l) = 0.0;
			Hvv(k, l) = 0.0;
			Svv(k, l) = 0.0;
			S(k, l) = 0.0;
			H(k, l) = 0.0;
		}
	}

	// find matrix elements for each term
	if (ka > 0)
		dpmAreaHessian2D(Ha, Sa);

	if (kl > 0)
		dpmPerimeterHessian2D(Hl, Sl);

	if (kb > 0)
		dpmBendingHessian2D(Hb,Sb);

	if (kc > 0)
		dpmRepulsiveHarmonicSprings2D(Hvv, Svv);

	// construct matrices
	for (k = 0; k < vertDOF; k++) {
		for (l = 0; l < vertDOF; l++) {
			H(k, l) = Ha(k, l) + Hl(k, l) + Hb(k, l) + Hvv(k, l);
			S(k, l) = -Sa(k, l) - Sl(k, l) - Sb(k, l) - Svv(k, l);
		}
	}
}

// construct hessian for area term
void dpm::dpmAreaHessian2D(Eigen::MatrixXd &Ha, Eigen::MatrixXd &Sa) {
	// local variables
	int nvtmp, ci, vim1, vi, vip1, vjm1, vj, vjp1;
	int kxm1, kx, kxp1, kym1, ky, kyp1, lxm1, lym1, lx, ly, lxp1, lyp1;
	double rho0, a0tmp, a02tmp, da, da_dxi, da_dyi, da_dxj, da_dyj;
	double lim1x, lix, liy, lim1y, ljm1x, ljm1y, ljx, ljy;

	// loop over cells
	rho0 = sqrt(a0[0]);
	for (ci = 0; ci < NCELLS; ci++) {
		// shape parameters for ci
		nvtmp = nv[ci];
		a0tmp = a0[ci];
		a02tmp = a0tmp * a0tmp;

		// fractional area strain
		da = (area(ci) / a0tmp) - 1.0;

		// loop over vertices
		for (vi = 0; vi < nvtmp; vi++) {
			// wrap vertices
			vim1 = (vi - 1 + nvtmp) % nvtmp;
			vip1 = (vi + 1) % nvtmp;

			// matrix indices
			kxm1 = NDIM * (gindex(ci, vim1));
			kym1 = NDIM * (gindex(ci, vim1)) + 1;

			kx = NDIM * (gindex(ci, vi));
			ky = NDIM * (gindex(ci, vi)) + 1;

			kxp1 = NDIM * (gindex(ci, vip1));
			kyp1 = NDIM * (gindex(ci, vip1)) + 1;

			// segment elements
			lim1x = x[kx] - x[kxm1];
			lim1y = x[ky] - x[kym1];

			lix = x[kxp1] - x[kx];
			liy = x[kyp1] - x[ky];

			if (pbc[0]) {
				lim1x -= L[0] * round(lim1x / L[0]);
				lix -= L[0] * round(lix / L[0]);
			}
			if (pbc[1]) {
				lim1y -= L[1] * round(lim1y / L[1]);
				liy -= L[1] * round(liy / L[1]);
			}

			// stress matrix
			Sa(kx, kyp1) = 0.5 * da * ((rho0 * rho0) / a0tmp);
			Sa(ky, kxp1) = -0.5 * da * ((rho0 * rho0) / a0tmp);

			Sa(kyp1, kx) = Sa(kx, kyp1);
			Sa(kxp1, ky) = Sa(ky, kxp1);

			// area derivatives (for stiffness matrix)
			da_dxi = 0.5 * (liy + lim1y);
			da_dyi = -0.5 * (lim1x + lix);

			// loop over other vertices, for area elasticity stiffness matrix
			for (vj = vi; vj < nvtmp; vj++) {

				// wrap jp1 and jm1
				vjp1 = (vj + 1) % nvtmp;
				vjm1 = (vj - 1 + nvtmp) % nvtmp;

				// dof elements
				lxm1 = NDIM * (gindex(ci, vjm1));
				lym1 = lxm1 + 1;

				lx = NDIM * (gindex(ci, vj));
				ly = lx + 1;

				lxp1 = NDIM * (gindex(ci, vjp1));
				lyp1 = lxp1 + 1;

				// j segments
				ljm1x = x[lx] - x[lxm1];
				if (pbc[0])
					ljm1x -= L[0] * round(ljm1x / L[0]);

				ljm1y = x[ly] - x[lym1];
				if (pbc[1])
					ljm1y -= L[1] * round(ljm1y / L[1]);

				ljx = x[lxp1] - x[lx];
				if (pbc[0])
					ljx -= L[0] * round(ljx / L[0]);

				ljy = x[lyp1] - x[ly];
				if (pbc[1])
					ljy -= L[1] * round(ljy / L[1]);

				// area derivatives
				da_dxj = 0.5 * (ljy + ljm1y);
				da_dyj = -0.5 * (ljm1x + ljx);

				// stiffness matrix
				Ha(kx, lx) = da_dxi * da_dxj * ((rho0 * rho0) / a02tmp);
				Ha(kx, ly) = da_dxi * da_dyj * ((rho0 * rho0) / a02tmp);

				Ha(ky, lx) = da_dyi * da_dxj * ((rho0 * rho0) / a02tmp);
				Ha(ky, ly) = da_dyi * da_dyj * ((rho0 * rho0) / a02tmp);

				Ha(lx, kx) = Ha(kx, lx);
				Ha(ly, kx) = Ha(kx, ly);

				Ha(lx, ky) = Ha(ky, lx);
				Ha(ly, ky) = Ha(ky, ly);
			}
		}
	}
}

// construct hessian for perimeter term
void dpm::dpmPerimeterHessian2D(Eigen::MatrixXd &Hl, Eigen::MatrixXd &Sl) {
	// local variables
	int nvtmp, ci, vim1, vi, vip1;
	int kxm1, kx, kxp1, kym1, ky, kyp1;
	double lim1x, lim1y, lix, liy, lim1, li, dlim1, dli, ulim1x, ulim1y, ulix, uliy;
	double l0im1, l0im1_sq, l0i, l0i_sq;
	double rho0, Kl;

	// loop over cells
	rho0 = sqrt(a0[0]);
	for (ci = 0; ci < NCELLS; ci++) {
		// number of vertices
		nvtmp = nv[ci];

		// prefactor scaled by length, will come out as dimensionless
		Kl = kl * (rho0 * rho0);

		for (vi = 0; vi < nvtmp; vi++) {
			// wrap vertices
			vim1 = (vi - 1 + nvtmp) % nvtmp;
			vip1 = (vi + 1) % nvtmp;

			// matrix indices
			kxm1 = NDIM * (gindex(ci, vim1));
			kym1 = NDIM * (gindex(ci, vim1)) + 1;

			kx = NDIM * (gindex(ci, vi));
			ky = NDIM * (gindex(ci, vi)) + 1;

			kxp1 = NDIM * (gindex(ci, vip1));
			kyp1 = NDIM * (gindex(ci, vip1)) + 1;

			// segment elements
			lim1x = x[kx] - x[kxm1];
			lim1y = x[ky] - x[kym1];

			lix = x[kxp1] - x[kx];
			liy = x[kyp1] - x[ky];

			if (pbc[0]) {
				lim1x -= L[0] * round(lim1x / L[0]);
				lix -= L[0] * round(lix / L[0]);
			}
			if (pbc[1]) {
				lim1y -= L[1] * round(lim1y / L[1]);
				liy -= L[1] * round(liy / L[1]);
			}

			// segment lengths
			lim1 = sqrt(lim1x * lim1x + lim1y * lim1y);
			li = sqrt(lix * lix + liy * liy);

			// segment strains
			l0im1 = l0[gindex(ci, vim1)];
			l0i = l0[gindex(ci, vi)];

			dlim1 = (lim1 / l0im1) - 1.0;
			dli = (li / l0i) - 1.0;

			l0im1_sq = l0im1 * l0im1;
			l0i_sq = l0i * l0i;

			// -- PERIMETER SPRINGS

			// unit vectors
			ulim1x = lim1x / lim1;
			ulim1y = lim1y / lim1;

			ulix = lix / li;
			uliy = liy / li;

			// 	STIFFNESS MATRIX

			// main diagonal
			Hl(kx, kx) = Kl * ((ulix * ulix) / l0i_sq + (ulim1x * ulim1x) / l0im1_sq);
			Hl(ky, ky) = Kl * ((uliy * uliy) / l0i_sq + (ulim1y * ulim1y) / l0im1_sq);

			Hl(kx, ky) = Kl * ((ulix * uliy) / l0i_sq + (ulim1x * ulim1y) / l0im1_sq);
			Hl(ky, kx) = Hl(kx, ky);

			// 1off diagonal
			Hl(kx, kxp1) = -Kl * (ulix * ulix) / l0i_sq;
			Hl(ky, kyp1) = -Kl * (uliy * uliy) / l0i_sq;

			Hl(kx, kyp1) = -Kl * (ulix * uliy) / l0i_sq;
			Hl(ky, kxp1) = Hl(kx, kyp1);

			// enforce symmetry in lower triangle
			Hl(kxp1, kx) = Hl(kx, kxp1);
			Hl(kyp1, ky) = Hl(ky, kyp1);

			Hl(kyp1, kx) = Hl(kx, kyp1);
			Hl(kxp1, ky) = Hl(ky, kxp1);

			// 	STRESS MATRIX

			// main diagonal
			Sl(kx, kx) = Kl * (dlim1 * ((ulim1y * ulim1y) / (l0im1 * lim1)) + dli * ((uliy * uliy) / (l0i * li)));
			Sl(ky, ky) = Kl * (dlim1 * ((ulim1x * ulim1x) / (l0im1 * lim1)) + dli * ((ulix * ulix) / (l0i * li)));

			Sl(kx, ky) = -Kl * (dlim1 * ((ulim1x * ulim1y) / (l0im1 * lim1)) + dli * ((ulix * uliy) / (l0i * li)));
			Sl(ky, kx) = Sl(kx, ky);

			// 1off diagonal
			Sl(kx, kxp1) = -Kl * dli * ((uliy * uliy) / (l0i * li));
			Sl(ky, kyp1) = -Kl * dli * ((ulix * ulix) / (l0i * li));

			Sl(kx, kyp1) = Kl * dli * ((ulix * uliy) / (l0i * li));
			Sl(ky, kxp1) = Sl(kx, kyp1);

			// enforce symmetry in lower triangle
			Sl(kxp1, kx) = Sl(kx, kxp1);
			Sl(kyp1, ky) = Sl(ky, kyp1);

			Sl(kyp1, kx) = Sl(kx, kyp1);
			Sl(kxp1, ky) = Sl(ky, kxp1);
		}
	}
}

// construct hessian for angular bending energy
void dpm::dpmBendingHessian2D(Eigen::MatrixXd &Hb, Eigen::MatrixXd &Sb){
	// local variables
	int gi, kxm1, kym1, kx, ky, kxp1, kyp1, kxp2, kyp2;
	double lxim1, lyim1, lx, ly, ulx, uly, ltmp, si, ci, th;
	double dtiim1x, dtiim1y, dtiix, dtiiy, dtiip1x, dtiip1y;
	double dtip1ix, dtip1iy, dtip1ip2x, dtip1ip2y, dtip1ip1x, dtip1ip1y;
	double dtim1, dti, dtip1;

	// non-dimensionalization
	double rho0, Kb;
	rho0 = sqrt(a0.at(0));
	Kb = kb * rho0 * rho0;

	// compute segment unit vectors matrix elements, segment lengths, normals & angle strain
	vector<double> ulxy(NVTOT,0.0);
	vector<double> uld(NVTOT,0.0);
	vector<double> l(NVTOT,0.0);
	vector<double> nx(NVTOT,0.0);
	vector<double> ny(NVTOT,0.0);
	vector<double> dth(NVTOT,0.0);
	for (gi=0; gi<NVTOT; gi++){
		// indexing
		kxm1 = NDIM*im1[gi];
		kym1 = kxm1 + 1;

		kx = NDIM*gi;
		ky = kx + 1;

		kxp1 = NDIM*ip1[gi];
		kyp1 = kxp1 + 1;

		// segment elements
		lxim1 = x[kx] - x[kxm1];
		lyim1 = x[ky] - x[kym1];

		lx = x[kxp1] - x[kx];
		ly = x[kyp1] - x[ky];

		// check periodic boundary conditions
		if (pbc[0]) {
			lxim1 -= L[0]*round(lxim1/L[0]);
			lx -= L[0]*round(lx/L[0]);
		}
		if (pbc[1]) {
			lyim1 -= L[1]*round(lyim1/L[1]);
			ly -= L[1]*round(ly/L[1]);
		}

		// segment length
		ltmp = sqrt(lx*lx + ly*ly);

		// unit vector elements
		ulx = lx/ltmp;
		uly = ly/ltmp;

		// unit vector matrix elements
		ulxy[gi] = ulx*uly;
		uld[gi] = 2.0*(uly*uly) - 1.0;

		// save length
		l[gi] = ltmp;

		// save normals
		nx[gi] = uly/ltmp;
		ny[gi] = -ulx/ltmp;

		// compute angles
		si = lx*lyim1 - ly*lxim1;
		ci = lx*lxim1 + ly*lyim1;
		th = atan2(si,ci);
		dth[gi] = th - t0[gi];
	}

	// loop over vertices
	rho0 = sqrt(a0[0]);
	for (gi=0; gi<NVTOT; gi++){
		// indexing
		kx = NDIM*gi;
		ky = kx + 1;

		kxp1 = NDIM*ip1[gi];
		kyp1 = kxp1 + 1;

		kxp2 = NDIM*ip1[ip1[gi]];
		kyp2 = kxp2 + 1;

		// first derivatives
	    dtiim1x = nx[im1[gi]];
	    dtiim1y = ny[im1[gi]];
	    
	    dtiip1x = nx[gi];
	    dtiip1y = ny[gi];
	    
	    dtiix = -(dtiim1x + dtiip1x);
	    dtiiy = -(dtiim1y + dtiip1y);

	    dtip1ix = dtiip1x;
	    dtip1iy = dtiip1y;

	    dtip1ip2x = nx[ip1[gi]];
	    dtip1ip2y = ny[ip1[gi]];

	    dtip1ip1x = -(dtip1ix + dtip1ip2x);
	    dtip1ip1y = -(dtip1iy + dtip1ip2y);

	    // angle strains
	    dtim1 = dth[im1[gi]];
	    dti = dth[gi];
	    dtip1 = dth[ip1[gi]];

		// -- Stiffness Matrix
    
	    // main diagonal block
	    Hb(kx,kx) = Kb * (pow(dtiim1x,2.0) + pow(dtiix,2.0) + pow(dtiip1x,2.0));
	    Hb(ky,ky) = Kb * (pow(dtiim1y,2.0) + pow(dtiiy,2.0) + pow(dtiip1y,2.0));
	    Hb(kx,ky) = Kb * ((dtiim1x*dtiim1y) + (dtiix*dtiiy) + (dtiip1x*dtiip1y));
	    Hb(ky,kx) = Hb(kx,ky);
	    
	    // 1off diagonal (enforce symmetry)
	    Hb(kx,kxp1) = Kb * (dtiix*dtiip1x + dtip1ix*dtip1ip1x);
	    Hb(ky,kyp1) = Kb * (dtiiy*dtiip1y + dtip1iy*dtip1ip1y);
	    Hb(kx,kyp1) = Kb * (dtiix*dtiip1y + dtip1ix*dtip1ip1y);
	    Hb(ky,kxp1) = Kb * (dtiiy*dtiip1x + dtip1iy*dtip1ip1x);
	    
	    Hb(kxp1,kx) = Hb(kx,kxp1);
	    Hb(kyp1,ky) = Hb(ky,kyp1);
	    Hb(kyp1,kx) = Hb(kx,kyp1);
	    Hb(kxp1,ky) = Hb(ky,kxp1);
	    
	    // 2off diagonal (enforce symmetry)
	    Hb(kx,kxp2) = Kb * dtip1ix * dtip1ip2x;
	    Hb(ky,kyp2) = Kb * dtip1iy * dtip1ip2y;
	    Hb(kx,kyp2) = Kb * dtip1ix * dtip1ip2y;
	    Hb(ky,kxp2) = Kb * dtip1iy * dtip1ip2x;
	    
	    Hb(kxp2,kx) = Hb(kx,kxp2);
	    Hb(kyp2,ky) = Hb(ky,kyp2);
	    Hb(kyp2,kx) = Hb(kx,kyp2);
	    Hb(kxp2,ky) = Hb(ky,kxp2);
	    
	    
	    // -- Stress Matrix
	        
	    // main diagonal block
	    Sb(kx,kx)       = (2.0*Kb*(dtip1 - dti)*ulxy[gi]/pow(l[gi],2.0)) + (2.0*Kb*(dti - dtim1)*ulxy[im1[gi]]/(pow(l[im1[gi]],2.0)));
	    Sb(ky,ky)       = -Sb(kx,kx);
	    Sb(kx,ky)       = Kb*(dtip1 - dti)*uld[gi]/pow(l[gi],2.0) + Kb*(dti - dtim1)*uld[im1[gi]]/pow(l[im1[gi]],2.0);
	    Sb(ky,kx)       = Sb(kx,ky);
	    
	    // 1off diagonal
	    Sb(kx,kxp1)     = 2.0*Kb*(dti - dtip1)*ulxy[gi]/pow(l[gi],2.0);
	    Sb(ky,kyp1)     = -Sb(kx,kxp1);
	    Sb(kx,kyp1)     = Kb*(dti - dtip1)*uld[gi]/pow(l[gi],2.0);
	    Sb(ky,kxp1)     = Sb(kx,kyp1);
	    
	    // enforce symmetry in lower triangle
	    Sb(kxp1,kx)     = Sb(kx,kxp1);
	    Sb(kyp1,ky)     = Sb(ky,kyp1);
	    Sb(kyp1,kx)     = Sb(kx,kyp1);
	    Sb(kxp1,ky)     = Sb(ky,kxp1);
	}

	// clear vectors from memory
	ulxy.clear();
	uld.clear();
	l.clear();
	nx.clear();
	ny.clear();
	dth.clear();
}

// construct hessian for interaction term
void dpm::dpmRepulsiveHarmonicSprings2D(Eigen::MatrixXd &Hvv, Eigen::MatrixXd &Svv) {
	// local variables
	int ci, cj, vi, vj, gi, gj;
	int mxi, myi, mxj, myj;
	double rho0, sij, dx, dy, rij, kij, h, uxij, uyij;

	// loop over cell pairs
	rho0 = sqrt(a0[0]);
	for (ci = 0; ci < NCELLS; ci++) {
		for (cj = ci + 1; cj < NCELLS; cj++) {

			// check if pair of cells is contact, only proceed if true
			if (cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2] > 0) {

				// loop over pairs of vertices on both cells, check for overlap, compute matrix elements
				for (vi = 0; vi < nv[ci]; vi++) {

					// matrix element indices (cell ci, vertex vi)
					gi = gindex(ci, vi);
					mxi = NDIM * gi;
					myi = mxi + 1;

					for (vj = 0; vj < nv[cj]; vj++) {
						// matrix element indices (cell cj, vertex vj)
						gj = gindex(cj, vj);
						mxj = NDIM * gj;
						myj = mxj + 1;

						// contact distance
						sij = r[gi] + r[gj];

						// get distance between vertices
						// particle distance
						dx = x[mxj] - x[mxi];
						if (pbc[0])
							dx -= L[0] * round(dx / L[0]);
						if (dx < sij) {
							dy = x[myj] - x[myi];
							if (pbc[1])
								dy -= L[1] * round(dy / L[1]);
							if (dy < sij) {
								rij = sqrt(dx * dx + dy * dy);
								if (rij < sij) {
									// spring constant
									kij = (kc * rho0 * rho0) / (sij * rij);

									// dimensionless overlap
									h = rij / sij;

									// derivatives of distance w.r.t. coordinates
									uxij = dx / rij;
									uyij = dy / rij;

									// compute stiffness and stress matrices (off diagonal, enforce symmetry in lower triangles)

									// -- stiffness matrix
									Hvv(mxi, mxj) = -((kc * rho0 * rho0) / (sij * sij)) * (uxij * uxij);
									Hvv(myi, myj) = -((kc * rho0 * rho0) / (sij * sij)) * (uyij * uyij);
									Hvv(mxi, myj) = -((kc * rho0 * rho0) / (sij * sij)) * (uxij * uyij);
									Hvv(myi, mxj) = -((kc * rho0 * rho0) / (sij * sij)) * (uyij * uxij);

									Hvv(mxj, mxi) = Hvv(mxi, mxj);
									Hvv(myj, myi) = Hvv(myi, myj);
									Hvv(mxj, myi) = Hvv(myi, mxj);
									Hvv(myj, mxi) = Hvv(mxi, myj);

									// -- stress matrix
									Svv(mxi, mxj) = kij * (1.0 - h) * (uyij * uyij);
									Svv(myi, myj) = kij * (1.0 - h) * (uxij * uxij);
									Svv(mxi, myj) = -kij * (1.0 - h) * (uxij * uyij);
									Svv(myi, mxj) = -kij * (1.0 - h) * (uxij * uyij);

									Svv(mxj, mxi) = Svv(mxi, mxj);
									Svv(myj, myi) = Svv(myi, myj);
									Svv(mxj, myi) = Svv(myi, mxj);
									Svv(myj, mxi) = Svv(mxi, myj);

									// add to diagonal, using off diagonals and reciprocity

									// -- stiffness matrix
									Hvv(mxi, mxi) -= Hvv(mxi, mxj);
									Hvv(myi, myi) -= Hvv(myi, myj);
									Hvv(mxi, myi) -= Hvv(mxi, myj);
									Hvv(myi, mxi) -= Hvv(myi, mxj);

									Hvv(mxj, mxj) -= Hvv(mxi, mxj);
									Hvv(myj, myj) -= Hvv(myi, myj);
									Hvv(mxj, myj) -= Hvv(mxi, myj);
									Hvv(myj, mxj) -= Hvv(myi, mxj);

									// -- stress matrix
									Svv(mxi, mxi) -= Svv(mxi, mxj);
									Svv(myi, myi) -= Svv(myi, myj);
									Svv(mxi, myi) -= Svv(mxi, myj);
									Svv(myi, mxi) -= Svv(myi, mxj);

									Svv(mxj, mxj) -= Svv(mxi, mxj);
									Svv(myj, myj) -= Svv(myi, myj);
									Svv(mxj, myj) -= Svv(mxi, myj);
									Svv(myj, mxj) -= Svv(myi, mxj);
								}
							}
						}
					}
				}
			}
		}
	}
}

/******************************

	P R I N T   T O

	C O N S O L E  &  F I L E

*******************************/

void dpm::printNeighborList() {
	int i,j,gi,d,boxid,sbtmp;
	double xtmp;

	cout << "NBX = " << NBX << endl;
	cout << "L_x = " << L[0] << ",  L_y = " << L[1] << endl;
	cout << "sb_x = " << sb[0] << ",  sb_y = " << sb[1] << endl;
	cout << "lb_x = " << lb[0] << ",  lb_y = " << lb[1] << endl;
	for (i=0; i<NBX; i++){
		cout << "nn[" << i << "]:  ";
		for (j=0; j<NNN; j++){
			cout << nn.at(i).at(j) << "  ";
		}
		cout << endl;
	}

	sortNeighborLinkedList2D();

	// loop over particles, plot location
	cout << endl << endl;
	cout << "PARTICLE BOX ID:" << endl;
	for (gi=0; gi<NVTOT; gi++){
		// 1. get cell id of current particle position
		boxid = 0;
		sbtmp = 1;
		for (d = 0; d < NDIM; d++) {
			// current location
			xtmp = x[NDIM * gi + d];

			// check out-of-bounds
			if (xtmp < 0){
				if (pbc[d])
					xtmp += L[d];
				else
					xtmp = 0.00001;
			}
			else if (xtmp > L[d]){
				if (pbc[d])
					xtmp -= L[d];
				else
					xtmp = 0.99999*L[d];
			}

			// add d index to 1d list
			boxid += floor(xtmp / lb[d]) * sbtmp;

			// increment dimensional factor
			sbtmp *= sb[d];
		}
		cout << "gi=" << gi << ";  boxid=" << boxid << "/" << NBX << endl;
		if (boxid < 0 || boxid >= NBX)
			cout << "** NOTE: boxid=" << boxid << ", BUT NBX=" << NBX << endl;
	}


	// loop over list
	cout << endl << endl;
	cout << "PARTICLE LIST ID:" << endl;
	for (gi=0; gi<NVTOT; gi++){
		cout << "gi=" << gi << ",  pi=" << gi+1 << ", list[pi]=" << list[gi+1] << endl;
	}
}

void dpm::printContactMatrix() {
	int ci, cj;

	for (ci = 0; ci < NCELLS; ci++) {
		for (cj = 0; cj < NCELLS; cj++) {
			if (ci > cj)
				cout << setw(5) << cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];
			else if (ci < cj)
				cout << setw(5) << cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
			else
				cout << setw(5) << 0;
		}
		cout << endl;
	}
}

void dpm::printConfiguration2D() {
	// local variables
	int ci, cj, vi, gi, ctmp, zc, zv;
	double xi, yi, dx, dy, Lx, Ly;

	// check if pos object is open
	if (!posout.is_open()) {
		cerr << "** ERROR: in printConfiguration2D, posout is not open, but function call will try to use. Ending here." << endl;
		exit(1);
	}
	else
		cout << "** In printConfiguration2D, printing particle positions to file..." << endl;

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

	// print stress info
	posout << setw(w) << left << "STRSS";
	posout << setw(wnum) << setprecision(pnum) << left << stress.at(0);
	posout << setw(wnum) << setprecision(pnum) << left << stress.at(1);
	posout << setw(wnum) << setprecision(pnum) << left << stress.at(2);
	posout << endl;

	// print coordinate for rest of the cells
	for (ci = 0; ci < NCELLS; ci++) {
		// get cell contact data
		zc = 0;
		zv = 0;
		for (cj = 0; cj < NCELLS; cj++) {
			if (ci != cj) {
				// contact info from entry ci, cj
				if (ci < cj)
					ctmp = cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2];
				else
					ctmp = cij[NCELLS * cj + ci - (cj + 1) * (cj + 2) / 2];

				// add to contact information
				zv += ctmp;
				if (ctmp > 0)
					zc++;
			}
		}

		// cell information
		posout << setw(w) << left << "CINFO";
		posout << setw(w) << left << nv.at(ci);
		posout << setw(w) << left << zc;
		posout << setw(w) << left << zv;
		posout << setw(wnum) << left << a0.at(ci);
		posout << setw(wnum) << left << area(ci);
		posout << setw(wnum) << left << perimeter(ci);
		posout << endl;

		// get initial vertex positions
		gi = gindex(ci, 0);
		xi = x.at(NDIM * gi);
		yi = x.at(NDIM * gi + 1);

		// place back in box center
		xi = fmod(xi, Lx);
		yi = fmod(yi, Ly);

		posout << setw(w) << left << "VINFO";
		posout << setw(w) << left << ci;
		posout << setw(w) << left << 0;

		// output initial vertex information
		posout << setw(wnum) << setprecision(pnum) << right << xi;
		posout << setw(wnum) << setprecision(pnum) << right << yi;
		posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
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
			posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
			posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR"
		   << " " << endl;
}

void dpm::printHessianEigenvalues2D(ofstream &hessout, Eigen::MatrixXd &M) {
	// check if pos object is open
	if (!hessout.is_open()) {
		cerr << "** ERROR: in printMatrixEigenvalues2D, hessout is not open, but function call will try to use. Ending here." << endl;
		exit(1);
	}
	else
		cout << "** In printMatrixEigenvalues2D, printing particle positions to file..." << endl;

	// compute eigenvalues from matrix, plot
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> dynamicalMatrixEigenmodes(M);

	// print to file
	hessout << vertDOF << endl;
	hessout << dynamicalMatrixEigenmodes.eigenvalues() << endl;
}
