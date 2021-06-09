/*

	BASIC FUNCTION DEFINITIONS for DPM class

	Jack Treado, 04/10/21

*/


#include "dpm.h"

// namespace
using namespace Eigen;
using namespace std;



/******************************

	C O N S T R U C T O R S  & 

		D E S T R U C T O R

*******************************/



// Main constructor
dpm::dpm(int n, int ndim, int seed){
	// local variables
	int d, i;

	// print to console
	cout << "** Instantiating configobj2D object, NCELLS = " << n << ",  ndim = " << ndim << ", seed = " << seed << " ..." << endl;

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

	// default boundary variables
	L.resize(NDIM);
	pbc.resize(NDIM);
	for (d=0; d<NDIM; d++){
		L[d] = 1.0;
		pbc[d] = 1;
	}

	// preferred area for each cell
	a0.resize(NCELLS);

	// macroscopic stress vector
	stress.resize(NDIM*(NDIM+1)/2);
	for (i=0; i<NDIM*(NDIM+1)/2; i++)
		stress.at(i) = 0.0;

	// contact network vector
	cij.resize(NCELLS*(NCELLS-1)/2);
	for (i=0; i<NCELLS*(NCELLS-1)/2; i++)
		cij.at(i) = 0; 

	// initialize nearest neighbor info
	NBX = 0;

	// seed random number generator
	srand48(seed);
}

// destructor
dpm::~dpm(){
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
	for (int i=0; i<NBX; i++)
		nn.at(i).clear();
	nn.clear();
	head.clear();
	last.clear();
	list.clear();

	
	if (posout.is_open())
		posout.close();

	if (hessout.is_open())
		hessout.close();

	if (ctcout.is_open())
		ctcout.close();
}



/******************************

	C E L L   S H A P E

	G E T T E R S

*******************************/


// get global vertex index gi given input cell index ci and vertex index vi
int dpm::gindex(int ci, int vi){
	return szList[ci] + vi;
} 


// get cell index ci and vertex index 
void dpm::cindices(int& ci, int& vi, int gi){
	for (int i=NCELLS-1; i>=0; i--){
		if (gi >= szList[i]){
			ci = i;
			vi = gi - szList[ci];
			break;
		}
	}
}


// get cell area
double dpm::area(int ci){
	// local variables
	int vi, vip1, gi, gip1, nvtmp;
	double dx, dy, xi, yi, xip1, yip1, areaVal = 0.0;

	// initial position: vi = 0
	nvtmp 	= nv.at(ci);
	gi 		= gindex(ci,0);
	xi 		= x[NDIM*gi];
	yi 		= x[NDIM*gi + 1];

	// loop over vertices of cell ci, get area by shoe-string method
	for (vi=0; vi<nvtmp; vi++){
		// next vertex
		gip1 = ip1[gi];
		gi++;

		// get positions (check minimum images)
		dx = x[NDIM*gip1] - xi;
		if (pbc[0])
			dx -= L[0]*round(dx/L[0]);
		xip1 = xi + dx;

		dy = x[NDIM*gip1 + 1] - yi;
		if (pbc[1])
			dy -= L[1]*round(dy/L[1]);
		yip1 = yi + dy;

		// increment area
		areaVal += xi*yip1 - xip1*yi;

		// set next coordinates
		xi = xip1;
		yi = yip1;
	}
	areaVal *= 0.5;

	return abs(areaVal);
}


// get cell perimeter
double dpm::perimeter(int ci){
	// local variables
	int vi, gi, gip1, nvtmp;
	double dx, dy, xi, yi, xip1, yip1, l, perimVal = 0.0;

	// initial position: vi = 0
	nvtmp 	= nv.at(ci);
	gi 		= gindex(ci,0);
	xi 		= x[NDIM*gi];
	yi 		= x[NDIM*gi + 1];

	// loop over vertices of cell ci, get perimeter
	for (vi=0; vi<nvtmp; vi++){
		// next vertex
		gip1 = ip1[gi];
		gi++;

		// get positions (check minimum images)
		dx = x[NDIM*gip1] - xi;
		if (pbc[0])
			dx -= L[0]*round(dx/L[0]);
		xip1 = xi + dx;

		dy = x[NDIM*gip1 + 1] - yi;
		if (pbc[1])
			dy -= L[1]*round(dy/L[1]);
		yip1 = yi + dy;

		// compute segment length
		l = sqrt(dx*dx + dy*dy);

		// add to perimeter
		perimVal += l;

		// update coordinates
		xi = xip1;
		yi = yip1;
	}

	// return perimeter
	return perimVal;
}


// get configuration packing fraction
double dpm::vertexPackingFraction2D(){
	int ci;
	double val, boxV, areaSum = 0.0;

	// numerator
	for (ci=0; ci<NCELLS; ci++)
		areaSum += area(ci) + 0.25*PI*pow(2.0*r.at(szList[ci]),2.0)*(0.5*nv.at(ci) - 1);

	// denominator
	boxV = L[0]*L[1];
	
	// return packing fraction
	val = areaSum/boxV;
	return val;
}

// get configuration "preferred" packing fraction
double dpm::vertexPreferredPackingFraction2D(){
	int ci;
	double val, boxV, areaSum = 0.0;

	// numerator
	for (ci=0; ci<NCELLS; ci++)
		areaSum += a0[ci] + 0.25*PI*pow(2.0*r.at(szList[ci]),2.0)*(0.5*nv.at(ci) - 1);

	// denominator
	boxV = L[0]*L[1];
	
	// return packing fraction
	val = areaSum/boxV;
	return val;
}




/******************************

	I N I T I A L -

			I Z A T I O N

*******************************/


// initialize vertex indexing
void dpm::initializeVertexIndexing2D(){
	int gi, vi, vip1, vim1, ci;

	// check that vertDOF has been assigned
	if (NVTOT <= 0){
		cout << "	** ERROR: in initializeVertexIndexing2D, NVTOT not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cout << "	** ERROR: in initializeVertexIndexing2D, vertDOF not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cout << "	** ERROR: in initializeVertexIndexing2D, nv vector not assigned. Need to initialize x, v, and F vectors in this function, so ending here." << endl;
		exit(1);
	}

	// save list of adjacent vertices
	im1.resize(NVTOT);
	ip1.resize(NVTOT);
	for (ci=0; ci<NCELLS; ci++){
		// vertex indexing
		for (vi=0; vi<nv.at(ci); vi++){
			// wrap local indices
			vim1 = (vi - 1 + nv.at(ci)) % nv.at(ci);
			vip1 = (vi + 1) % nv.at(ci);

			// get global wrapped indices
			gi 			= gindex(ci,vi);
			im1.at(gi) 	= gindex(ci,vim1);
			ip1.at(gi) 	= gindex(ci,vip1);
		}
	}

	// initialize vertex configuration vectors
	x.resize(vertDOF);
	v.resize(vertDOF);
	F.resize(vertDOF);
}


// initialize vertex shape parameters based on nv
void dpm::initializeVertexShapeParameters(double calA0, int nref){
	// local variables
	int gi, ci, vi, nvtmp;
	double rtmp, calA0tmp, calAntmp;

	// check that vertDOF has been assigned
	if (NVTOT <= 0){
		cout << "	** ERROR: in initializeVertexShapeParameters, NVTOT not assigned. Ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cout << "	** ERROR: in initializeVertexShapeParameters, vertDOF not assigned. Ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cout << "	** ERROR: in initializeVertexShapeParameters, nv vector not assigned. Ending here." << endl;
		exit(1);
	}

	// resize shape paramters
	l0.resize(NVTOT);
	t0.resize(NVTOT);
	r.resize(NVTOT);

	// loop over cells, determine shape parameters
	for (ci=0; ci<NCELLS; ci++){
		// number of vertices on cell ci
		nvtmp 		= nv.at(ci);

		// a0 based on nv
		rtmp 		= (double)nvtmp/nref;
		a0.at(ci) 	= rtmp*rtmp;

		// shape parameter
		calAntmp 	= nvtmp*tan(PI/nvtmp)/PI;
		calA0tmp 	= calA0*calAntmp;

		// l0 and vertex radii
		gi 			= szList.at(ci);
		for (vi=0; vi<nv.at(ci); vi++){
			l0.at(gi+vi) 	= 2.0*sqrt(PI*calA0tmp*a0.at(ci))/nvtmp;
			t0.at(gi+vi) 	= 0.0;
			r.at(gi+vi) 	= 0.5*l0.at(gi+vi);
		}
	}
}


// initialize bidisperse cell system, single calA0
void dpm::bidisperse2D(double calA0, int nsmall, double smallfrac, double sizefrac){
	// local variables
	double calA0tmp, calAntmp, rtmp, areaSum;
	int vim1, vip1, gi, ci, vi, nlarge, smallN, largeN, NVSMALL;

	// print to console
	cout << "** initializing bidisperse DPM particles in 2D ..." << endl;

	// number of vertices on large particles
	nlarge = round(sizefrac*nsmall);

	// total number of vertices
	smallN 	= round(smallfrac*NCELLS);
	largeN 	= NCELLS - smallN;
	NVSMALL = nsmall*smallN;
	NVTOT 	= NVSMALL + nlarge*largeN;
	vertDOF = NDIM*NVTOT;

	// szList and nv (keep track of global vertex indices)
	nv.resize(NCELLS);
	szList.resize(NCELLS);

	nv.at(0) = nsmall;
	for (ci=1; ci<NCELLS; ci++){
		if (ci < smallN){
			nv.at(ci) = nsmall;
			szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);
		}
		else{
			nv.at(ci) = nlarge;
			szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);
		}
	}

	// initialize vertex shape parameters
	initializeVertexShapeParameters(calA0, nsmall);

	// initialize vertex indexing
	initializeVertexIndexing2D();
}

// initialize gaussian polydisperse cell system, single calA0
void dpm::gaussian2D(double dispersion, double calA0, int n1){
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
	for (ci=1; ci<NCELLS; ci++){
		// use Box-Muller to generate polydisperse sample
		r1 = drand48();
		r2 = drand48();
		grv = sqrt(-2.0*log(r1))*cos(2.0*PI*r2);
		nvtmp = floor(dispersion*n1*grv + n1);
		if (nvtmp < nvmin)
			nvtmp = nvmin;

		// store size of cell ci
		nv.at(ci) = nvtmp;
		szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

		// add to total NV count
		NVTOT += nvtmp;
	}
	vertDOF = NDIM*NVTOT;

	// initialize vertex shape parameters
	initializeVertexShapeParameters(calA0, n1);

	// initialize vertex indexing
	initializeVertexIndexing2D();
}


// set sinusoidal preferred angle
void dpm::sinusoidalPreferredAngle(double thA, double thK){
	int ci, vi, gi;
	double thR;

	// print to console
	cout << "** setting initial th0 values to sinusoids, thA = " << thA << ", thK = " << thK << " ..." << endl;

	// loop over cells
	gi = 0;
	for (ci=0; ci<NCELLS; ci++){
		thR = (2.0*PI)/nv.at(ci);
		for (vi=0; vi<nv.at(ci); vi++){
			t0.at(gi) = thA*thR*sin(thR*thK*vi);
			gi++;
		}
	}
}

// initialize CoM positions of cells using SP FIRE
void dpm::initializePositions2D(double phi0, double Ftol){
	// local variables
	int i, d, ci, cj, vi, vj, gi, cellDOF = NDIM*NCELLS;
	double areaSum, xtra = 1.1;

	// local disk vectors
	vector<double> drad(NCELLS,0.0);
	vector<double> dpos(cellDOF,0.0);
	vector<double> dv(cellDOF,0.0);
	vector<double> dF(cellDOF,0.0);

	// print to console
	cout << "** initializing particle positions using 2D SP model and FIRE relaxation ..." << endl;

	// initialize box size based on packing fraction
	areaSum = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		areaSum += a0.at(ci) + 0.25*PI*pow(l0.at(ci),2.0)*(0.5*nv.at(ci) - 1);

	// set box size
	for (d=0; d<NDIM; d++)
		L.at(d) = pow(areaSum/phi0,1.0/NDIM);

	// initialize cell centers randomly
	for (ci=0; ci<cellDOF; ci += 2)
		dpos.at(ci) = L[ci % 2]*drand48();
	for (ci=cellDOF-1; ci>0; ci -= 2)
		dpos.at(ci) = L[ci % 2]*drand48();

	// set radii of SP disks
	for (ci=0; ci<NCELLS; ci++)
		drad.at(ci) = xtra*sqrt((2.0*a0.at(ci))/(nv.at(ci)*sin(2.0*PI/nv.at(ci))));

	// FIRE VARIABLES
	double P  		= 0;	
	double fnorm 	= 0;
	double vnorm 	= 0;
	double alpha   	= alpha0;

	double dt0 		= 1e-2;
	double dtmax   	= 10*dt0;
	double dtmin   	= 1e-8*dt0;

	int npPos      	= 0;
	int npNeg      	= 0;

	int fireit    	= 0;
	int itmax 		= 1e6;
	double fcheck  	= 10*Ftol;

	// interaction variables
	double rij, sij, dtmp, ftmp, vftmp;
	double dr[NDIM];

	// initial step size
	dt = dt0;


	// loop until force relaxes
	while ((fcheck > Ftol) && fireit < itmax){
		// FIRE step 1. Compute P
		P = 0.0;
		for (i=0; i<cellDOF; i++)
			P += dv[i]*dF[i];

		// FIRE step 2. adjust simulation based on net motion of degrees of freedom
		if (P > 0){
			// increase positive counter
			npPos++;

			// reset negative counter
			npNeg = 0;

			// alter simulation if enough positive steps have been taken
			if (npPos > NMIN){
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
				cout << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
				exit(1);
			}

			// take half step backwards, reset velocities
			for (i=0; i<cellDOF; i++){
				// take half step backwards
				dpos[i] -= 0.5*dt*dv[i];

				// reset velocities
				dv[i] = 0.0;
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

		// FIRE step 3. First VV update
		for (i=0; i<cellDOF; i++)
			dv[i] += 0.5*dt*dF[i];

		// FIRE step 4. adjust velocity magnitude
		fnorm = 0.0;
		vnorm = 0.0;
		for (i=0; i<cellDOF; i++){
			fnorm 	+= dF[i]*dF[i];
			vnorm 	+= dv[i]*dv[i];
		}
		fnorm = sqrt(fnorm);
		vnorm = sqrt(vnorm);
		if (fnorm > 0){
			for (i=0; i<cellDOF; i++)
				dv[i] = (1 - alpha)*dv[i] + alpha*(vnorm/fnorm)*dF[i];
		}

		// FIRE step 4. Second VV update
		for (i=0; i<cellDOF; i++){
			dpos[i] += dt*dv[i];
			dF[i] = 0.0;
		}

		// FIRE step 5. Update forces
		for (ci=0; ci<NCELLS; ci++){
			for (cj=ci+1; cj<NCELLS; cj++){

				// contact distance
				sij = drad[ci] + drad[cj];

				// true distance
				rij = 0.0;
				for (d=0; d<NDIM; d++){
					// get distance element
					dtmp = dpos[NDIM*cj + d] - dpos[NDIM*ci + d];
					if (pbc[d])
						dtmp -= L[d]*round(dtmp/L[d]);

					// add to true distance
					rij += dtmp*dtmp;

					// save in distance array
					dr[d] = dtmp;
				}
				rij = sqrt(rij);

				// check distances
				if (rij < sij){
					// force magnitude
					ftmp = kc*(1.0 - (rij/sij))/sij;
					
					// add to vectorial force
					for (d=0; d<NDIM; d++){
						vftmp = ftmp*(dr[d]/rij);
						dF[NDIM*ci + d] -= vftmp;
						dF[NDIM*cj + d] += vftmp;
					}
				}
			}
		}

		// FIRE step 5. Final VV update
		for (i=0; i<cellDOF; i++)
			dv[i] += 0.5*dt*dF[i];

		// update forces to check
		fcheck = 0.0;
		for (i=0; i<cellDOF; i++)
			fcheck += dF[i]*dF[i];
		fcheck = sqrt(fcheck/NCELLS);

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
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
			cout << "	** Pdir = " << P/(fnorm*vnorm) << endl;
			cout << "	** alpha = " << alpha << endl;
		}

		// update iterate
		fireit++;
	}
	// check if FIRE converged
	if (fireit == itmax){
		cout << "	** FIRE minimization did not converge, fireit = " << fireit << ", itmax = " << itmax << "; ending." << endl;
		exit(1);
	}
	else{
		cout << endl << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===========================================" << endl;
		cout << " 	F I R E 						" << endl;
		cout << "		M I N I M I Z A T I O N 	" << endl;
		cout << "	C O N V E R G E D! 				" << endl << endl;

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
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<nv.at(ci); vi++){
			// get global vertex index
			gi = gindex(ci,vi);

			// length from center to vertex
			dtmp = sqrt((2.0*a0.at(ci))/(nv.at(ci)*sin((2.0*PI)/nv.at(ci))));

			// set positions
			x.at(NDIM*gi) 		= dtmp*cos((2.0*PI*vi)/nv.at(ci)) + dpos.at(NDIM*ci) + 1e-2*l0[gi]*drand48();
			x.at(NDIM*gi + 1)	= dtmp*sin((2.0*PI*vi)/nv.at(ci)) + dpos.at(NDIM*ci + 1) + 1e-2*l0[gi]*drand48();
		}
	}
}


// initialize neighbor linked list
void dpm::initializeNeighborLinkedList2D(double boxLengthScale){
	// local variables
	double llscale = 2.0*r.at(0);
	int i, d, nntmp, scx;

	// print to console
	cout << "** initializing neighbor linked list, boxLengthScale = " << boxLengthScale;

	// initialize box length vectors
	NBX = 1;
	sb.resize(NDIM);
	lb.resize(NDIM);
	for (d=0; d<NDIM; d++){
		// determine number of cells along given dimension by rmax
		sb[d] = round(L[d]/(boxLengthScale*llscale));

		// just in case, if < 3, change to 3 so box neighbor checking will work
		if (sb[d] < 3)
			sb[d] = 3;

		// determine box length by number of cells
		lb[d] = L[d]/sb[d];

		// count total number of cells
		NBX *= sb[d];
	}

	// initialize list of box nearest neighbors
	scx = sb[0];
	nn.resize(NBX);

	// loop over cells, save forward neighbors for each box
	for (i=0; i<NBX; i++){
		// reshape entry
		nn[i].resize(NNN);
		
		// neighbors
		nn[i][0] 			= (i + 1) % NBX; 			// right neighbor (i+1)
		nn[i][1] 			= (i + scx) % NBX;			// top neighbor (j+1)
		nntmp 				= (i + NBX - scx) % NBX;	// bottom neighbor (j-1)
		nn[i][2] 			= (nn[i][1] + 1) % NBX;		// top-right neighbor (i+1, j+1)
		nn[i][3] 			= nntmp + 1;				// bottom-right neighbor (i+1, j-1)

		// right-hand bc (periodic)
		if ((i+1) % scx == 0){
			nn[i][0] = i - scx + 1;
			nn[i][2] = nn[i][1]  - scx + 1;
			nn[i][3] = nntmp - scx + 1;
		}
	}

	// linked-list variables
	head.resize(NBX);
	last.resize(NBX);
	list.resize(NVTOT+1);

	// print box info to console
	cout << ";  initially NBX = " << NBX << " ..." << endl;
}




/******************************

	E D I T I N G   &

			U P D A T I N G

*******************************/


// sort vertices into neighbor linked list
void dpm::sortNeighborLinkedList2D(){
	// local variables
	int d, gi, boxid, sbtmp;

	// reset linked list info
	fill(list.begin(), list.end(), 0);
	fill(head.begin(), head.end(), 0);
	fill(last.begin(), last.end(), 0);

	// sort vertices into linked list
	for (gi=0; gi<NVTOT; gi++){
		// 1. get cell id of current particle position
		boxid = 0;
		sbtmp = 1;
		for (d=0; d<NDIM; d++){
			// add d index to 1d list
			boxid += floor(x[NDIM*gi + d]/lb[d])*sbtmp;

			// increment dimensional factor
			sbtmp *= sb[d];
		}

		// 2. add to head list or link within list
		// NOTE: particle ids are labelled starting from 1, setting to 0 means end of linked list
		if (head[boxid] == 0){
			head[boxid] = gi + 1;
			last[boxid] = gi + 1;
		}
		else{
			list[last[boxid]] = gi + 1;
			last[boxid] = gi + 1;
		}
	}
}


// change size of particles
void dpm::scaleParticleSizes2D(double scaleFactor){
	// local variables
	int gi, ci, vi, xind, yind;
	double xi, yi, cx, cy, dx, dy;

	// loop over cells, scale
 	for (ci=0; ci<NCELLS; ci++){
		// scale preferred area
		a0[ci] *= scaleFactor*scaleFactor;

		// first global index for ci
		gi = szList.at(ci);

		// compute cell center of mass
		xi = x[NDIM*gi];
		yi = x[NDIM*gi + 1];
		cx = xi; 
		cy = yi;
		for (vi=1; vi<nv.at(ci); vi++){
			dx = x.at(NDIM*(gi+vi)) - xi;
			if (pbc[0])
				dx -= L[0]*round(dx/L[0]);

			dy = x.at(NDIM*(gi+vi) + 1) - yi;
			if (pbc[1])
				dy -= L[1]*round(dy/L[1]);

			xi += dx;
			yi += dy;

			cx += xi;
			cy += yi;
		}
		cx /= nv.at(ci);
		cy /= nv.at(ci);

		for (vi=0; vi<nv.at(ci); vi++){
			// x and y inds
			xind = NDIM*(gi+vi);
			yind = xind + 1;

			// closest relative position
			dx = x[xind] - cx;
			if (pbc[0])
				dx -= L[0]*round(dx/L[0]);

			dy = x[yind] - cy;
			if (pbc[1])
				dy -= L[1]*round(dy/L[1]);

			// update vertex positions
			x[xind] 		+= (scaleFactor - 1.0)*dx;
			x[yind] 		+= (scaleFactor - 1.0)*dy;

			// scale vertex radii
			r[gi+vi] 		*= scaleFactor;
			l0[gi+vi] 		*= scaleFactor;
		}
	}
}


// remove rattlers from contact network, return number of rattlers
int dpm::removeRattlers(){
	// local variables
	int ci, cj, ctmp, rvv, rcc, nr, nm=1;

	// loop over rows, eliminate contacts to rattlers until nm = 0
	while (nm > 0){
		// reset number of rattlers
		nr = 0;

		// number of "marginal" rattlers to be removed
		nm = 0;
		for (ci=0; ci<NCELLS; ci++) {
			// get number of contacts on cell ci
			rvv = 0;
			rcc = 0;
			for (cj=0; cj<NCELLS; cj++){
				if (ci != cj){
					if (ci > cj)
						ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2];
					else
						ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2];
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

					for (cj=0; cj<NCELLS; cj++) {
						// delete contact between ci and cj
						if (ci != cj){
							if (ci > cj)
								cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2] = 0;
							else
								cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2] = 0; 
						}
					}
				}
			}
		}
	}

	// return total number of rattlers
	return nr;
}





/******************************

	D P M  F O R C E 

			U P D A T E S

*******************************/


void dpm::resetForcesAndEnergy(){
	fill(F.begin(), F.end(), 0.0);
	fill(stress.begin(), stress.end(), 0.0);
	U = 0.0;
}


void dpm::shapeForces2D(){
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
	for (gi=0; gi<NVTOT; gi++){

		// -- Area force (and get cell index ci)
		if (ci < NCELLS){
			if (gi == szList[ci]){
				// shape information
				nvtmp = nv[ci];
				a0tmp = a0[ci];

				// preferred segment length of last segment
				l0im1 = l0[im1[gi]];

				// compute area deviation
				atmp = area(ci);
				da = (atmp/a0tmp) - 1.0;

				// update potential energy
				U += 0.5*ka*(da*da);

				// shape force parameters
				fa = da*(rho0/a0tmp);
				fb = kb/rho0;

				// compute cell center of mass
				xi = x[NDIM*gi];
				yi = x[NDIM*gi + 1];
				cx = xi; 
				cy = yi;
				for (vi=1; vi<nvtmp; vi++){
					// get distances between vim1 and vi
					dx = x[NDIM*(gi+vi)] - xi;
					dy = x[NDIM*(gi+vi) + 1] - yi;
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// add to centers
					xi += dx;
					yi += dy;

					cx += xi;
					cy += yi;
				}
				cx /= nvtmp;
				cy /= nvtmp;

				// get coordinates relative to center of mass
				rix = x[NDIM*gi] - cx;
				riy = x[NDIM*gi + 1] - cy;

				// get prior adjacent vertices
				rim2x = x[NDIM*im1[im1[gi]]] - cx;
				rim2y = x[NDIM*im1[im1[gi]] + 1] - cy;
				if (pbc[0])
					rim2x -= L[0]*round(rim2x/L[0]);
				if (pbc[1])
					rim2y -= L[1]*round(rim2y/L[1]);

				rim1x = x[NDIM*im1[gi]] - cx;
				rim1y = x[NDIM*im1[gi] + 1] - cy;
				if (pbc[0])
					rim1x -= L[0]*round(rim1x/L[0]);
				if (pbc[1])
					rim1y -= L[1]*round(rim1y/L[1]);

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
		rip1x = x[NDIM*ip1[gi]] - cx;
		rip1y = x[NDIM*ip1[gi] + 1] - cy;
		if (pbc[0])
			rip1x -= L[0]*round(rip1x/L[0]);
		if (pbc[1])
			rip1y -= L[1]*round(rip1y/L[1]);


		// -- Area force
		F[NDIM*gi] 		+= 0.5*fa*(rim1y - rip1y);
		F[NDIM*gi + 1] 	+= 0.5*fa*(rip1x - rim1x);


		// -- Perimeter force
		lix 	= rip1x - rix;
		liy 	= rip1y - riy;

		// segment lengths
		lim1 	= sqrt(lim1x*lim1x + lim1y*lim1y);
		li 		= sqrt(lix*lix + liy*liy);

		// segment deviations
		dlim1  	= (lim1/l0im1) - 1.0;
		dli 	= (li/l0i) - 1.0;

		// segment forces
		flim1 	= kl*(rho0/l0im1);
		fli 	= kl*(rho0/l0i);

		// add to forces
		F[NDIM*gi] 		+= (fli*dli*lix/li) - (flim1*dlim1*lim1x/lim1);
		F[NDIM*gi + 1] 	+= (fli*dli*liy/li) - (flim1*dlim1*lim1y/lim1);
		
		// update potential energy
		U += 0.5*kl*(dli*dli);


		// -- Bending force
		if (kb > 0){
			// get ip2 for third angle
			rip2x = x[NDIM*ip1[ip1[gi]]] - cx;
			rip2y = x[NDIM*ip1[ip1[gi]] + 1] - cy;
			if (pbc[0])
				rip2x -= L[0]*round(rip2x/L[0]);
			if (pbc[1])
				rip2y -= L[1]*round(rip2y/L[1]);

			// get last segment length
			lip1x = rip2x - rip1x;
			lip1y = rip2y - rip1y;

			// get angles
			sinim1 = lim1x*lim2y - lim1y*lim2x;
			cosim1 = lim1x*lim2x + lim1y*lim2y;

			sini = lix*lim1y - liy*lim1x;
			cosi = lix*lim1x + liy*lim1y;

			sinip1 = lip1x*liy - lip1y*lix;
			cosip1 = lip1x*lix + lip1y*liy;

			// get normal vectors
			nim1x = lim1y;
			nim1y = -lim1x;

			nix = liy;
			niy = -lix;

			// get change in angles
			dtim1 = atan2(sinim1,cosim1) - t0[im1[gi]];
			dti = atan2(sini,cosi) - t0[gi];
			dtip1 = atan2(sinip1,cosip1) - t0[ip1[gi]];

			// get delta delta theta's
			ddtim1 = (dti - dtim1)/(lim1*lim1);
			ddti = (dti - dtip1)/(li*li);

			// add to force
			F[NDIM*gi] 		+= fb*(ddtim1*nim1x + ddti*nix);
			F[NDIM*gi + 1] 	+= fb*(ddtim1*nim1y + ddti*niy);

			// update potential energy
			U += 0.5*kb*(dti*dti);
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


void dpm::repulsiveVertexForces2D(){
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, dx, dy, rho0;
	double ftmp, fx, fy;

	// sort particles
	sortNeighborLinkedList2D();

	// get fundamental length
	rho0 = sqrt(a0[0]);

	// reset contact network
	fill(cij.begin(), cij.end(), 0);

	// loop over boxes in neighbor linked list
	for (bi=0; bi<NBX; bi++){

		// get start of list of vertices
		pi = head[bi];

		// loop over linked list
		while (pi > 0){
			// real particle index
			gi = pi - 1;

			// next particle in list
			pj = list[pi];

			// loop down neighbors of pi in same cell
			while (pj > 0){
				// real index of pj
				gj = pj - 1;

				if (gj == ip1[gi] || gj == im1[gi]){
					pj = list[pj];
					continue;
				}

				// contact distance
				sij = r[gi] + r[gj];

				// particle distance
				dx = x[NDIM*gj] - x[NDIM*gi];
				if (pbc[0])
					dx -= L[0]*round(dx/L[0]);
				if (dx < sij){
					dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);
					if (dy < sij){
						rij = sqrt(dx*dx + dy*dy);
						if (rij < sij){
							// force scale
							ftmp 				= kc*(1 - (rij/sij))*(rho0/sij);
							fx 					= ftmp*(dx/rij);
							fy 					= ftmp*(dy/rij);

							// add to forces
							F[NDIM*gi] 			-= fx;
							F[NDIM*gi + 1] 		-= fy;

							F[NDIM*gj] 			+= fx;
							F[NDIM*gj + 1] 		+= fy;

							// increae potential energy
							U += 0.5*kc*pow((1 - (rij/sij)),2.0);

							// add to virial stress
							stress[1] += dx*fx;
							stress[2] += dy*fy;
							stress[3] += 0.5*(dx*fy + dy*fx);

							// add to contacts
							cindices(ci, vi, gi);
							cindices(cj, vj, gj);

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
			for (bj=0; bj<NNN; bj++){
				// get first particle in neighboring cell
				pj = head[nn[bi][bj]];

				// loop down neighbors of pi in same cell
				while (pj > 0){
					// real index of pj
					gj = pj - 1;

					if (gj == ip1[gi] || gj == im1[gi]){
						pj = list[pj];
						continue;
					}
					// contact distance
					sij = r[gi] + r[gj];

					// particle distance
					dx = x[NDIM*gj] - x[NDIM*gi];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);
					if (dx < sij){
						dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
						if (pbc[1])
							dy -= L[1]*round(dy/L[1]);
						if (dy < sij){
							rij = sqrt(dx*dx + dy*dy);
							if (rij < sij){
								// force scale
								ftmp 				= kc*(1 - (rij/sij))*(rho0/sij);
								fx 					= ftmp*(dx/rij);
								fy 					= ftmp*(dy/rij);

								// add to forces
								F[NDIM*gi] 			-= fx;
								F[NDIM*gi + 1] 		-= fy;

								F[NDIM*gj] 			+= fx;
								F[NDIM*gj + 1] 		+= fy;

								// increae potential energy
								U += 0.5*kc*pow((1 - (rij/sij)),2.0);

								// add to virial stress
								stress[1] += dx*fx;
								stress[2] += dy*fy;
								stress[3] += 0.5*(dx*fy + dy*fx);

								// add to contacts
								cindices(ci, vi, gi);
								cindices(cj, vj, gj);

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
	stress[0] *= (rho0/(L[0]*L[1]));
	stress[1] *= (rho0/(L[0]*L[1]));
	stress[2] *= (rho0/(L[0]*L[1]));
}














/******************************

	D P M  

		I N T E G R A T O R S

*******************************/


void dpm::setdt(double dt0){
	// local variables
	int i;
	double ta, tl, tb, tmin, rho0;

	// typical length
	rho0 = sqrt(a0.at(0));

	// set typical time scales
	ta = rho0/sqrt(ka);
	tl = (rho0*l0.at(0))/sqrt(ka*kl);
	tb = (rho0*l0.at(0))/sqrt(ka*kb);

	// set main time scale as min
	tmin = 1e8;
	if (ta < tmin)
		tmin = ta;
	if (tl < tmin)
		tmin = tl;
	if (tb < tmin)
		tmin = tb;

	// set dt
	dt = dt0*tmin;
}


void dpm::vertexFIRE2D(double Ftol, double dt0){
	// local variables
	int i;
	double rho0;

	// FIRE variables
	double P, fnorm, fcheck, vnorm, alpha, dtmax, dtmin;
	int npPos, npNeg, fireit;

	// set dt based on geometric parameters
	setdt(dt0);

	// Initialize FIRE variables
	P  			= 0;	
	fnorm 		= 0;
	vnorm 		= 0;
	alpha   	= alpha0;

	dtmax   	= 10.0*dt;
	dtmin   	= 1e-2*dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;

	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);

	// length scale
	rho0 = sqrt(a0.at(0));

	// relax forces using FIRE
	while (fcheck > Ftol && fireit < itmax){
		// compute P
		P = 0.0;
		for (i=0; i<vertDOF; i++)
			P += v[i]*F[i];

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
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
		}

		// Adjust simulation based on net motion of degrees of freedom
		if (P > 0){
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
				cout << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
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

		// VV VELOCITY UPDATE #1
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*dt*F[i];

		// compute fnorm, vnorm and P
		fnorm = 0.0;
		vnorm = 0.0;
		for (i=0; i<vertDOF; i++){
			fnorm 	+= F[i]*F[i];
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
		for (i=0; i<vertDOF; i++){
			// update position
			x[i] += dt*v[i];

			// recenter in box
			if (x[i] > L[i % NDIM] && pbc[i % NDIM])
				x[i] -= L[i % NDIM];
			else if (x[i] < 0 && pbc[i % NDIM])
				x[i] += L[i % NDIM];
		}

		// update forces
		resetForcesAndEnergy();
		repulsiveVertexForces2D();
		shapeForces2D();

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*F[i]*dt;

		// update fcheck based on fnorm (= force per degree of freedom)
		fcheck = 0.0;
		for (i=0; i<vertDOF; i++)
			fcheck += F[i]*F[i];
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
		cout << endl << endl;
	}
}












/******************************

	D P M  

		P R O T O C O L S

*******************************/


void dpm::vertexCompress2Target2D(double Ftol, double dt0, double phi0Target, double dphi0){
	// local variables
	int it = 0, itmax = 1e4;
	double phi0 = vertexPreferredPackingFraction2D();
	double scaleFactor, P, Sxy;

	// loop while phi0 < phi0Target
	while (phi0 < phi0Target && it < itmax){
		// scale particle sizes
		scaleParticleSizes2D(scaleFactor);

		// update phi0
		phi0 = vertexPreferredPackingFraction2D();

		// relax configuration
		vertexFIRE2D(Ftol, dt0);

		// get scale factor
		scaleFactor = sqrt((phi0 + dphi0)/phi0);

		// get updated pressure
		P = 0.5*(stress[0] + stress[1]);
		Sxy = stress[2];

		// print to console
		cout << endl << endl;
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
		cout << endl << endl;

		// update iterate
		it++;
	}
}


void dpm::vertexJamming2D(double Ftol, double Ptol, double dt0, double dphi0, bool plotCompression){
	// local variables
	int k=0, nr;
	bool jammed, overcompressed, undercompressed;
	double pcheck, phi0, rH, r0, rL, rho0, scaleFactor;

	// initialize binary root search parameters
	r0 = sqrt(a0.at(0));
	rH = -1;
	rL = -1;

	// initialize preferred packing fraction
	phi0 = vertexPreferredPackingFraction2D();
	
	// save initial state
	vector<double> xsave(vertDOF,0.0);
	vector<double> rsave(vertDOF,0.0);
	vector<double> l0save(vertDOF,0.0);
	vector<double> t0save(vertDOF,0.0);
	vector<double> a0save(vertDOF,0.0);

	xsave = x;
	rsave = r;
	l0save = l0;
	t0save = t0;
	a0save = a0;


	// loop until jamming is found
	while (!jammed && k < itmax){
		// set length scale by 1st particle preferred area
		rho0 = sqrt(a0.at(0));

		// relax configuration
		vertexFIRE2D(Ftol, dt0);

		// update pressure
		pcheck = 0.5*(stress[0] + stress[1]);

		// remove rattlers
		nr = removeRattlers();

		// boolean checks for jamming
		undercompressed = ((pcheck < 2.0*Ptol && rH < 0) || (pcheck < Ptol && rH > 0));
		overcompressed = (pcheck > 2.0*Ptol);
		jammed = (pcheck < 2.0*Ptol && pcheck > Ptol && rH > 0 && rL > 0);

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << " 	Q U A S I S T A T I C  				" << endl;
		cout << " 	  	I S O T R O P I C 				" << endl;
		cout << "			C O M P R E S S I O N 		" << endl << endl;
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
		cout << "	* # of rattlers = " << nr << endl << endl;
		cout << "	* undercompressed = " << undercompressed << endl;
		cout << "	* overcompressed = " << overcompressed << endl;
		cout << "	* jammed = " << jammed << endl << endl;
		if (plotCompression)
			printConfiguration2D();
		printContactMatrix();
		cout << endl << endl;

		// update particle scaleFactor based on target check
		if (rH < 0){
			// if still undercompressed, then grow until overcompressed found
			if (undercompressed){
				r0 = rho0;
				scaleFactor = sqrt((phi0 + dphi0)/phi0);
			}
			// if first overcompressed, decompress by dphi/2 until unjamming
			else if (overcompressed){
				// current = upper bound length scale r
	            rH = rho0;

	            // save first overcompressed state
				r0 = rH;
				xsave = x;
				rsave = r;
				l0save = l0;
				t0save = t0;
				a0save = a0;

	            // shrink particle sizes
	            scaleFactor = sqrt((phi0 - 0.5*dphi0)/phi0);

	            // print to console
				cout << "	-- -- overcompressed for the first time, scaleFactor = " << scaleFactor << endl;
			}
		}
		else{
			if (rL < 0){
				// if first undercompressed, save last overcompressed state, begin root search
				if (undercompressed){
					// current = new lower bound length scale r
					rL = rho0;

					// load state
					x = xsave;
					r = rsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute new scale factor by root search
		            scaleFactor = 0.5*(rH + rL)/r0;

					// print to console
					cout << "	-- -- undercompressed for the first time, scaleFactor = " << scaleFactor << endl;
					cout << "	-- -- BEGINNING ROOT SEARCH IN ENTHALPY MIN PROTOCOL..." << endl;
				}
				// if still overcompressed, decrement again
				else if (overcompressed){
					// current = upper bound length scale r
		            rH = rho0;

		            // save overcompressed state
					r0 = rH;
					xsave = x;
					rsave = r;
					l0save = l0;
					t0save = t0;
					a0save = a0;

		            // keep shrinking at same rate until unjamming
		            scaleFactor = sqrt((phi0 - 0.5*dphi0)/phi0);

		            // print to console
					cout << "	-- -- overcompressed, still no unjamming, scaleFactor = " << scaleFactor << endl;
				}
			}
			else{
				// if found undercompressed state, go to state between undercompressed and last overcompressed states (from saved state)
				if (undercompressed){
					// current = new lower bound length scale r
					rL = rho0;

					// load state
					x = xsave;
					r = rsave;
					l0 = l0save;
					t0 = t0save;
					a0 = a0save;

					// compute new scale factor by root search
		            scaleFactor = 0.5*(rH + rL)/r0;

					// print to console
					cout << "	-- -- undercompressed, scaleFactor = " << scaleFactor << endl;

				}
				else if (overcompressed){
					// current = upper bound length scale r
		            rH = rho0;

		            // save overcompressed state
					r0 = rH;
					xsave = x;
					rsave = r;
					l0save = l0;
					t0save = t0;
					a0save = a0;

					// compute new scale factor
		            scaleFactor = 0.5*(rH + rL)/r0;

					// print to console
					cout << "	-- -- overcompressed, scaleFactor = " << scaleFactor << endl;
				}
				else if (jammed){
					cout << "	** At k = " << k << ", target pressure found!" << endl;
					cout << " WRITING ENTHALPY-MINIMIZED CONFIG TO FILE" << endl;
					cout << " ENDING COMPRESSION SIMULATION" << endl;
					scaleFactor = 1.0;
					if (!plotCompression)
						printConfiguration2D();
					break;
				}
			}
		}

		// scale particle sizes
		scaleParticleSizes2D(scaleFactor);

		// update packing fraction
		phi0 = vertexPreferredPackingFraction2D();
	}
}










/******************************

	P R I N T   T O

	C O N S O L E  &  F I L E

*******************************/

void dpm::printContactMatrix(){
	int ci, cj;

	for (ci=0; ci<NCELLS; ci++){
		for (cj=0; cj<NCELLS; cj++){
			if (ci > cj)
				cout << cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2] << "  ";
			else if (ci < cj)
				cout << cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2] << "  "; 
			else
				cout << "0  ";
		}
		cout << endl;
	}
}

void dpm::printConfiguration2D(){
	// local variables
	int ci, cj, vi, gi, ctmp, zc, zv;
	double xi, yi, dx, dy, Lx, Ly;

	// print to console
	cout << "** Printing particle positions to file" << endl;

	// save box sizes
	Lx = L.at(0);
	Ly = L.at(1);

	// print information starting information
	posout << setw(w) << left << "NEWFR" << " " << endl;
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
	for (ci=0; ci<NCELLS; ci++){

		// get cell contact data
		zc = 0;
		zv = 0;
		for (cj=0; cj<NCELLS; cj++){
			if (ci != cj){
				// contact info from entry ci, cj
				if (ci < cj)
					ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]; 
				else
					ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]; 

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
		gi = gindex(ci,0);
		xi = x.at(NDIM*gi);
		yi = x.at(NDIM*gi + 1);

		// place back in box center
		xi = fmod(xi,Lx);
		yi = fmod(yi,Ly);

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
		for (vi=1; vi<nv.at(ci); vi++){
			// get global vertex index for next vertex
			gi++;

			// get next vertex positions
			dx = x.at(NDIM*gi) - xi;
			if (pbc[0])
				dx -= Lx*round(dx/Lx);
			xi += dx;

			dy = x.at(NDIM*gi + 1) - yi;
			if (pbc[1])
				dy -= Ly*round(dy/Ly);
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
	posout << setw(w) << left << "ENDFR" << " " << endl;
}





