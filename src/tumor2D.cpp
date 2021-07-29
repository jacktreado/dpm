/*

	FUNCTION DEFINITIONS for tumor2D class

	Jack Treado, 06/04/21

	** TO DO 06/24/21
	1. Add adipocyte / tumor boundary sim initialization function
	2. Make protocol that compresses ONLY so that we have premade initial conditions, don't have to spend overhead
	3. Add constant pressure via piston attached to left-side wall, box vol. changes as cells invade
		* Need to add kinetic term to stress tensor

*/


#include "tumor2D.h"
#include "dpm.h"

// namespace
using namespace Eigen;
using namespace std;



/*********************************

	C O N S T R U C T O R S 

	&

	D E S T R U C T O R S 

**********************************/


// read-in constructor
tumor2D::tumor2D(string &inputFileStr,int seed) : dpm(2) {
	// open file
	ifstream inputobj(inputFileStr.c_str());
	if (!inputobj.is_open()){
		cerr << "** ERROR: In tumor2D constructor, could not open file " << inputFileStr << ", ending here. " << endl;
		exit(1);
	}

	// set variables to default
	gamtt=0.0; v0=0.0; Dr0=0.0; Ds=0.0; kecm = 0.0; ecmbreak = 0.0; pbc[0]=0; pbc[1]=1;

	// local variables
	int nvtmp, ci, vi, i;
	double val;
	double lxtmp, lytmp;
	double s1, s2, s3;
	double wp1, wp2;
	double a0tmp, l0tmp, t0tmp;
	double xtmp, ytmp, rtmp;
	string inputStr;

	// LINE 1: should be NEWFR
	getline(inputobj, inputStr);
	if (inputStr.compare(0,5,"NEWFR") != 0){
		cerr << "** ERROR: In tumor2D constructor, first line of input file NOT NEWFR, first line = " << inputStr << ". Ending." << endl;
		exit(1);
	}

	// read in simulation information from header
	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"NUMCL %d %d",&NCELLS,&tN);
	cout << "\t ** " << inputStr << endl;

	// verify input file
	if (NCELLS< 1){
		cerr << "** ERROR: in tumor2D constructor, NCELLStmp = " << NCELLS << ". Ending here." << endl;
		exit(1);
	}
	else if (tN < 1){
		cerr << "** ERROR: in tumor2D constructor, tNtmp = " << tN << ". Ending here." << endl;
		exit(1);
	}

	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"TSTEP %lf",&val);
	cout << "\t ** " << inputStr << endl;

	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"TCURR %lf",&val);
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

	// initialize stress
	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"STRSS %lf %lf %lf",&s1,&s2,&s3);
	cout << "\t ** " << inputStr << endl;

	stress.at(0) = s1;
	stress.at(1) = s2;
	stress.at(2) = s3;

	// initialize wall pressure
	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"WPRSS %lf %lf",&wp1,&wp2);
	cout << "\t ** " << inputStr << endl;

	wpress.resize(2);
	wpress.at(0) = wp1;
	wpress.at(1) = wp2;


	// szList and nv (keep track of global vertex indices)
	nv.resize(NCELLS);
	szList.resize(NCELLS);
	a0.resize(NCELLS);

	// initialize ecm + crawling variables
	psi.resize(tN);
	Dr.resize(tN);

	fill(psi.begin(), psi.end(), 0.0);
	fill(Dr.begin(), Dr.end(), 0.0);

	pinpos.resize(NDIM * (NCELLS - tN));
	pinattach.resize(NCELLS - tN);

	// initialize NVTOT to 0
	NVTOT = 0;

	// loop over cells, read in coordinates
	cout << "\n\n** LOOPING OVER DATA, PRINTING INFO..." << endl;
	for (ci=0; ci<NCELLS; ci++){
		// first parse cell info
		getline(inputobj, inputStr);
		if (ci < tN)
			sscanf(inputStr.c_str(),"CINFO %d %*d %*d %lf %*lf %*lf %*lf %*f",&nvtmp,&a0tmp);
		else
			sscanf(inputStr.c_str(),"CINFO %d %*d %*d %lf %*lf %*lf %*lf %*lf %*lf",&nvtmp,&a0tmp);

		// print to console
		cout << "\t ** " << inputStr << endl;

		// store in vectors
		nv.at(ci) = nvtmp;
		a0.at(ci) = a0tmp;

		// update NVTOT
		NVTOT += nvtmp;

		// increment szList
		if (ci > 0)
			szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

		// loop over vertices, store coordinates
		for (vi=0; vi<nvtmp; vi++){
			// parse vertex coordinate info
			getline(inputobj, inputStr);
			sscanf(inputStr.c_str(),"VINFO %*d %*d %lf %lf %lf %lf %lf",&xtmp,&ytmp,&rtmp,&l0tmp,&t0tmp);

			// print to console
			cout << "\t ** " << inputStr << endl;

			// push back 
			x.push_back(xtmp);
			x.push_back(ytmp);
			r.push_back(rtmp);
			l0.push_back(l0tmp);
			t0.push_back(t0tmp);
		}
	}
	vertDOF = NDIM * NVTOT;
	cout << "** FINISHED LOOPING OVER DATA, CLOSING INPUT FILE OBJ\n" << endl;

	// initialize contact network (NEED TO DO HERE, NOT DONE IN DPM CONSTRUCTOR)
	cij.resize(NCELLS * (NCELLS - 1) / 2);
	for (i = 0; i < NCELLS * (NCELLS - 1) / 2; i++)
		cij.at(i) = 0;

	// initialize vertex indexing
	initializeVertexIndexing2D();

	// set initial l0
	setl0_init();

	// close input file object
	inputobj.close();

	// seed random number generator
	srand48(seed);
}




/*********************************

	T U M O R   C E L L

	I N I T I A L I Z A T I O N

**********************************/


// initialize single tumor cell for growth as monolayer or single crawler
void tumor2D::initializeSingleTumorCell(){
	// local variables
	int vi;
	int dtmp;

	// only proceed if tN has not been set yet, but nv etc has
	if (tN != 0){
		cout << "\t ** ERROR: in initializeSingleTumorCell, tN = " << tN << " which is not 0, cannot reset. Ending here." << endl;
		exit(1);
	}
	if (NVTOT <= 0){
		cout << "	** ERROR: in initializeSingleTumorCell, NVTOT not assigned. Ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cout << "	** ERROR: in initializeSingleTumorCell, vertDOF not assigned. Ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cout << "	** ERROR: in initializeSingleTumorCell, nv vector not assigned. Ending here." << endl;
		exit(1);
	}

	// initialize tN to 1
	tN = 1;

	// create first cell centered at origin
	for (vi=0; vi<nv.at(0); vi++){
		// length from center to vertex
		dtmp = sqrt((2.0*a0.at(0))/(nv.at(0)*sin((2.0*PI)/nv.at(0))));

		// set positions
		x.at(NDIM*vi) 		= dtmp*cos((2.0*PI*vi)/nv.at(0)) + 1e-2*l0[vi]*drand48();
		x.at(NDIM*vi + 1)	= dtmp*sin((2.0*PI*vi)/nv.at(0)) + 1e-2*l0[vi]*drand48();
	}
}


// set l0_init to be l0
void tumor2D::setl0_init(){
	// check that l0 set
	if (NVTOT <= 0){
		cerr << "\t ** ERROR: in setl0_init, NVTOT = " << NVTOT << ", which is <= 0. Ending. " << endl;
		exit(1);
	}
	if (l0.size() != NVTOT){
		cerr << "\t ** ERROR: in setl0_init, l0 size = " << l0.size() << ", which is != " << NVTOT << ". Ending. " << endl;
		exit(1);
	}

	// resize l0_init
	l0_init.resize(NVTOT);

	// fill with l0
	for (int gi=0; gi<NVTOT; gi++)
		l0_init.at(gi) = l0.at(gi);
}


// set monolayer positions to center of box
void tumor2D::initializeTumorMonolayerPositions(double phi0, double Ftol, double kwell){
	// local variables
	int i, d, ci, cj, vi, vj, gi, cellDOF = NDIM * NCELLS;
	double areaSum, dnorm, xtra = 1.1;

	// initialize ecm + crawling variables
	psi.resize(tN);
	Dr.resize(tN);
	setl0_init();

	fill(psi.begin(), psi.end(), 2.0*PI*drand48());
	fill(Dr.begin(), Dr.end(), Dr0);

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

			// also add well force toward center of box
			dnorm = sqrt(pow(dpos[NDIM*ci] - 0.5*L[0],2.0) + pow(dpos[NDIM*ci + 1] - 0.5*L[1],2.0));
			dF[NDIM * ci] -= kwell*(dpos[NDIM*ci] - 0.5*L[0])/dnorm;
			dF[NDIM * ci + 1] -= kwell*(dpos[NDIM*ci + 1] - 0.5*L[1])/dnorm;
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
		cout << "	C O N V E R G E D! 				" << endl<< endl;

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


// Initialize collection of tumor cells and adipocytes
// 
// NOTE: parameters for adipocytes(tumors) 
// start with a(t) in camelCase variables
// 
// parameters
// -- calA0: 		shape parameter
// -- disp: 		size dispersity
// -- areaRatio: 	ratio of adipocyte preferred areas to tumor preferred areas
// -- NV: 			number of vertices on particle with average sqrt(a0)
// 
void tumor2D::initializeTumorInterface(double aCalA0, double tCalA0, double aDisp, double tDisp, double areaRatio, int aNV, int tNV){
	// local variables
	int ci, nvtmp;
	double lenscale, r1, r2, grv;

	// print to console
	cout << "** initializing gaussian tumor and adipocyte DPM particles in 2D with size dispersions:" << endl;
	cout << "\t** aDisp = " << aDisp << endl;
	cout << "\t** tDisp = " << tDisp << endl;
	cout << "** setting up nv + szList, setting shape parameters and initializing indexing ..." << endl;

	// szList and nv (keep track of global vertex indices)
	nv.resize(NCELLS);
	szList.resize(NCELLS);

	// initialize ecm + crawling variables
	psi.resize(tN);
	Dr.resize(tN);

	fill(psi.begin(), psi.end(), 0.0);
	fill(Dr.begin(), Dr.end(), 0.0);

	pinpos.resize(NDIM * (NCELLS - tN));
	pinattach.resize(NCELLS - tN);

	// initialize number of vertices on each cell
	nv.at(0) = tNV;
	NVTOT = tNV;
	for (ci=1; ci<NCELLS; ci++){
		// use Box-Muller to generate polydisperse sample
		r1 = drand48();
		r2 = drand48();
		grv = sqrt(-2.0*log(r1))*cos(2.0*PI*r2);
		if (ci < tN)
			nvtmp = floor(tDisp*tNV*grv + tNV);
		else
			nvtmp = floor(aDisp*aNV*grv + aNV);

		if (nvtmp < nvmin)
			nvtmp = nvmin;

		// store size of cell ci
		nv.at(ci) = nvtmp;
		szList.at(ci) = szList.at(ci-1) + nv.at(ci-1);

		// add to total NV count
		NVTOT += nvtmp;
	}
	vertDOF = NDIM * NVTOT;

	// resize shape paramters
	l0.resize(NVTOT);
	l0_init.resize(NVTOT);
	t0.resize(NVTOT);
	r.resize(NVTOT);

	// initialize particle sizes based on areaRatio
	for (ci=0; ci<NCELLS; ci++){
		if (ci < tN){
			lenscale = (double) nv.at(ci) / tNV;
			initializeVertexShapeParameters(ci,tCalA0,lenscale);
		}
		else{
			lenscale = (double) (nv.at(ci) * sqrt(areaRatio)) / aNV;
			initializeVertexShapeParameters(ci,aCalA0,lenscale);
		}
	}

	// initialize l0_init
	setl0_init();

	// initialize vertex indexing
	initializeVertexIndexing2D();
}


void tumor2D::initializeTumorInterfacePositions(double phi0, double Ftol, double prt){
	// local variables
	int i, d, ci, cj, vi, vj, gi, cellDOF = NDIM * NCELLS;
	double areaSum, xtra = 1.05, xi, Ldiv;

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
	L.at(1) = sqrt(areaSum/(2.0*phi0));
	L.at(0) = 2.0*L[1];

	// dividing wall position between adipocytes and 
	Ldiv = prt * L[0];

	// initialize tumor cell centers in left-hand partition of the box
	for (ci=0; ci<tN; ci++){
		dpos.at(NDIM*ci) 		= (Ldiv - 2.0*drad[ci])*drand48() + drad[ci];
		dpos.at(NDIM*ci + 1) 	= L[1]*drand48();
	}

	// initialize WAT cell centers to the right
	for (ci=tN; ci<NCELLS; ci++){
		dpos.at(NDIM*ci) 		= (L[0] - Ldiv - 2.0*drad[ci])*drand48() + Ldiv + drad[ci];
		dpos.at(NDIM*ci + 1) 	= L[1]*drand48();
	}

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

			// x boundary forces
			xi = dpos[NDIM * ci];
			if (ci < tN) {
				if (xi < drad[ci])
					dF[NDIM*ci] += (1.0 - (xi/drad[ci]))/drad[ci];
				else if (xi > Ldiv - drad[ci])
					dF[NDIM*ci] -= (1.0 - ((Ldiv - xi)/drad[ci]))/drad[ci];
			}
			else {
				if (xi < Ldiv + drad[ci])
					dF[NDIM*ci] += (1.0 - ((xi - Ldiv)/drad[ci]))/drad[ci];
				else if (xi > L[0] - drad[ci])
					dF[NDIM*ci] -= (1.0 - ((L[0] - xi)/drad[ci]))/drad[ci];
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


/*********************************

	E D I T I N G   &

			U P D A T I N G

**********************************/

// divide a single cell, assume preallocation
void tumor2D::divide(int ci){
	// local variable
	int gi, g1tmp, g2tmp, vi, ci1, ci2, nv1, nv2, icut1, icut2, nh1, nh2, vitmp;
	double Dx, Dy, Dnorm, dhatx, dhaty, nhatx, nhaty, xtmp, ytmp;

	// cell indices
	ci1 = tN - 1;
	ci2 = ci1 + 1;

	// number of vertices on each cell
	nv1 	= nv.at(ci1);
	nv2 	= nv.at(ci2);

	// end point of stitch
	icut1 	= szList.at(ci1) + floor(0.5*nv1) - 1;
	icut2 	= szList.at(ci1) + nv1 - 1;

	// get size of semicircles
	nh1 	= floor(0.5*nv1);
	nh2 	= nv1 - nh1;

	// move second-half of vertices from 1 -> first half of 2
	for (gi=0; gi<nh2; gi++){
		// get temporary global vertex indices
		g1tmp = icut1 + gi + 1;
		g2tmp = szList.at(ci2) + gi;

		// make old vertices part of cell 2 
		x.at(NDIM*g2tmp) = x.at(NDIM*g1tmp);
		x.at(NDIM*g2tmp + 1) = x.at(NDIM*g1tmp + 1);
	}

	// get centerline
	Dx = x.at(NDIM*szList.at(ci1)) - x.at(NDIM*icut1);
	if (pbc[0])
		Dx -= L[0]*round(Dx/L[0]);
	Dy = x.at(NDIM*szList.at(ci1) + 1) - x.at(NDIM*icut1 + 1);
	if (pbc[1])
		Dy -= L[1]*round(Dy/L[1]);
	Dnorm = sqrt(Dx*Dx + Dy*Dy);

	// centerline unit vector
	dhatx = Dx/Dnorm;
	dhaty = Dy/Dnorm;

	// centerline normal vector
	nhatx = -dhaty;
	nhaty = dhatx;

	// stitch vertices on centerline (cell 1)
	vitmp = 1;
	// for (gi=nh1+1; gi<nv1; gi++){
	// 	// get new position
	// 	xtmp += 
	// }
}







/******************************

	B I O L O G I C A L  

	F U N C T I O N S

*******************************/


// -- CRAWLING CELLS

// update psi based on Dr
void tumor2D::psiDiffusion(){
	// local variables
	int ci;
	double r1, r2, grv;

	// update director for each cell
	for (ci=0; ci<tN; ci++){
		// generate random variable
		r1 = drand48();
		r2 = drand48();
		grv = sqrt(-2.0*log(r1))*cos(2.0*PI*r2);

		// update director for cell ci
		psi[ci] += sqrt(2.0*dt*Dr[ci])*grv;
	}
}


void tumor2D::activeBrownianCrawlerUpdate(){
	// local variables
	int gi, ci, vi;
	double cx, cy, rix, riy, ux, uy, psitmp, dpsi, v0tmp, rnorm;

	// loop over all cells, vertices
	gi = 0; 
	for (ci=0; ci<tN; ci++){
		// center of mass of cell
		com2D(ci,cx,cy);

		// loop over vertices
		for (vi=0; vi<nv[ci]; vi++){
			// x position
			rix = x[NDIM*gi] - cx;
			if (pbc[0])
				rix -= L[0]*round(rix/L[0]);

			// y position
			riy = x[NDIM*gi + 1] - cy;
			if (pbc[1])
				riy -= L[1]*round(riy/L[1]);

			// get angular distance from psi
			psitmp = atan2(riy,rix);
			dpsi = psitmp - psi[ci];
			dpsi -= 2.0*PI*round(dpsi/(2.0*PI));

			// get velocity scale
			v0tmp = v0*(0.01 + 0.99*exp(-pow(dpsi,2.0)/(2.0*Ds*Ds)));

			// get unit vectors
			rnorm = sqrt(rix*rix + riy*riy);
			ux = rix/rnorm;
			uy = riy/rnorm;

			// add to forces
			F[NDIM*gi] += v0tmp*ux;
			F[NDIM*gi + 1] += v0tmp*uy;

			// update global vertex
			gi++;
		}
	}
}



// -- ADIPOCYTE ECM ATTACHEMENT

// update 
void tumor2D::updateECMAttachments(bool attach){
	// check that vectors have been created
	if (pinattach.size() != NCELLS-tN){
		cerr << "	** ERROR: in updatePinAttachments, number of pins != number of adipocytes, so ending here." << endl;
		exit(1);
	}
	if (pinpos.size() != NDIM*(NCELLS - tN)){
		cerr << "	** ERROR: in updatePinAttachments, number of pins != number of adipocytes, so ending here." << endl;
		exit(1);
	}

	// local variables
	int ci;
	double cx, cy;

	// set all attached
	fill(pinattach.begin(), pinattach.end(), attach);

	// update positions
	for (ci=tN; ci<NCELLS; ci++){
		com2D(ci,cx,cy);
		pinpos.at(NDIM*(ci-tN)) = cx;
		pinpos.at(NDIM*(ci-tN) + 1) = cy;
	}
}

// add to force
void tumor2D::adipocyteECMAdhesionForces(){
	int gi, ci, vi, xind, yind, nvtmp, pinCellInd;
	double dpinx, dpiny, dpin, cx, cy;

	// loop over cells
	for (ci=tN; ci<NCELLS; ci++){
		// get pin indices
		pinCellInd = ci - tN;

		// tmp number of vertices
		nvtmp = nv[ci];

		// get distance to pin if attached
		if (pinattach[pinCellInd]){
			// get centers of mass
			com2D(ci,cx,cy);

			// get distance to pin
			dpinx = pinpos[NDIM*pinCellInd] - cx;
			if (pbc[0])
				dpinx -= L[0]*round(dpinx/L[0]);

			dpiny = pinpos[NDIM*pinCellInd + 1] - cy;
			if (pbc[1])
				dpiny -= L[1]*round(dpiny/L[1]);

			dpin = sqrt(pow(dpinx,2.0) + pow(dpiny,2.0));

			// pin forces on each vertex
			gi = gindex(ci,0);
			for (vi=0; vi<nvtmp; vi++){
				// positions using global indexing
				xind = NDIM*(gi+vi);
				yind = xind + 1;

				// if an adipocyte and pin is intact, compute force due to pinning spring
				if (dpin < ecmbreak*sqrt(a0[ci])){
					F[xind] 	+= (kecm/nvtmp)*dpinx;
					F[yind] 	+= (kecm/nvtmp)*dpiny;
				}
				else
					pinattach[pinCellInd] = 0;
			}
		}
	}
}



/******************************

	F O R C E 

	U P D A T E S

*******************************/


void tumor2D::resetForcesAndEnergy(){
	fill(F.begin(), F.end(), 0.0);
	fill(stress.begin(), stress.end(), 0.0);
	fill(wpress.begin(), wpress.end(), 0.0);
	U = 0.0;
}



void tumor2D::repulsiveTumorForces() {
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, xij, dx, dy, rho0, xi, yi, ri;
	double ftmp, fx, fy;

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
							// scaled distance
							xij = rij/sij;

							// force magnitude
							ftmp = kc*(1 - xij)/sij;

							// increase potential energy
							U += 0.5*kc*pow(1.0 - xij,2.0);

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
								// scaled distance
								xij = rij/sij;

								// force magnitude
								ftmp = kc*(1 - xij)/sij;

								// increase potential energy
								U += 0.5*kc*pow(1.0 - xij,2.0);

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

void tumor2D::stickyTumorForces() {
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, dx, dy, rho0, xi, yi, ri;
	double ftmp, fx, fy;
	vector<int> ztt(NVTOT,0);

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
							// only tumor cells attract
							if (rij > cutij && ci < tN && cj < tN){
								// force scale
								ftmp = kint*(xij - 1.0 - l2)/sij;

								// increase potential energy
								U += -0.5*kint*pow(1.0 + l2 - xij,2.0);
							}
							else{
								if ((ci < tN && cj < tN) || rij < sij){
									// force scale
									ftmp = kc*(1 - xij)/sij;

									// increase potential energy
									U += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
								}
								else{
									pj = list[pj];
									continue;
								}
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

							// if both tumor cells, add to contact list for surface tension
							if (ci < tN && cj < tN){
								ztt[gi]++;
								ztt[gj]++;
							}
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
								// only tumor cells attract
								if (rij > cutij && ci < tN && cj < tN){
									// force scale
									ftmp = kint*(xij - 1.0 - l2)/sij;

									// increase potential energy
									U += -0.5*kint*pow(1.0 + l2 - xij,2.0);
								}
								else{
									if ((ci < tN && cj < tN) || xij < 1.0){
										// force scale
										ftmp = kc*(1 - xij)/sij;

										// increase potential energy
										U += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
									}
									else{
										pj = list[pj];
										continue;
									}
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

								// if both tumor cells, add to contact list for surface tension
								if (ci < tN && cj < tN){
									ztt[gi]++;
									ztt[gj]++;
								}
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

	// update l0 based on ztt 
	for (gi=0; gi<NVTOT; gi++){
		if (ztt[gi] > 0 && ztt[ip1[gi]] > 0)
			l0[gi] = l0_init[gi]*(1.0 - (gamtt/kl));
		else
			l0[gi] = l0_init[gi];
	}
}

void tumor2D::repulsiveTumorInterfaceForces() {
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, xij, dx, dy, rho0, xi, yi, ri;
	double ftmp, fx, fy;

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

			// check boundary forces
			xi = x[NDIM*gi];
			ri = r[gi];
			if (xi < ri){
				// update forces
				fx = kc*(1.0 - (xi/ri))/ri;
				F[NDIM*gi] += fx;

				// update wall stress
				wpress[0] += fx/L[1];
			}
			else if (xi > L[0] - ri){
				// update forces
				fx = -kc*(1.0 - ((L[0] - xi)/ri))/ri;
				F[NDIM*gi] += fx;

				// update wall stresses
				wpress[1] -= fx/L[1];
			}

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
							// scaled distance
							xij = rij/sij;

							// force magnitude
							ftmp = kc*(1 - xij)/sij;

							// increase potential energy
							U += 0.5*kc*pow(1.0 - xij,2.0);

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
								// scaled distance
								xij = rij/sij;

								// force magnitude
								ftmp = kc*(1 - xij)/sij;

								// increase potential energy
								U += 0.5*kc*pow(1.0 - xij,2.0);

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

void tumor2D::stickyTumorInterfaceForces(){
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, dx, dy, rho0, xi, yi, ri;
	double ftmp, fx, fy;
	vector<int> ztt(NVTOT,0);

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
							// only tumor cells attract
							if (rij > cutij && ci < tN && cj < tN){
								// force scale
								ftmp = kint*(xij - 1.0 - l2)/sij;

								// increase potential energy
								U += -0.5*kint*pow(1.0 + l2 - xij,2.0);
							}
							else{
								if ((ci < tN && cj < tN) || rij < sij){
									// force scale
									ftmp = kc*(1 - xij)/sij;

									// increase potential energy
									U += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
								}
								else{
									pj = list[pj];
									continue;
								}
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

							// if both tumor cells, add to contact list for surface tension
							if (ci < tN && cj < tN){
								ztt[gi]++;
								ztt[gj]++;
							}
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
								// only tumor cells attract
								if (rij > cutij && ci < tN && cj < tN){
									// force scale
									ftmp = kint*(xij - 1.0 - l2)/sij;

									// increase potential energy
									U += -0.5*kint*pow(1.0 + l2 - xij,2.0);
								}
								else{
									if ((ci < tN && cj < tN) || xij < 1.0){
										// force scale
										ftmp = kc*(1 - xij)/sij;

										// increase potential energy
										U += 0.5*kc*(pow(1.0 - xij,2.0) - l1*l2);
									}
									else{
										pj = list[pj];
										continue;
									}
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

								// if both tumor cells, add to contact list for surface tension
								if (ci < tN && cj < tN){
									ztt[gi]++;
									ztt[gj]++;
								}
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

	// update l0 based on ztt 
	for (gi=0; gi<NVTOT; gi++){
		if (ztt[gi] == 0 && ztt[ip1[gi]] == 0)
			l0[gi] = l0_init[gi]*(1.0 - (gamtt/kl));
		else
			l0[gi] = l0_init[gi];
	}
}






void tumor2D::repulsiveTumorForceUpdate() {
	resetForcesAndEnergy();
	repulsiveTumorForces();
	shapeForces2D();
}

void tumor2D::stickyTumorForceUpdate() {
	resetForcesAndEnergy();
	stickyTumorForces();
	shapeForces2D();
}

void tumor2D::repulsiveTumorInterfaceForceUpdate() {
	resetForcesAndEnergy();
	repulsiveTumorInterfaceForces();
	shapeForces2D();
}

void tumor2D::stickyTumorInterfaceForceUpdate() {
	resetForcesAndEnergy();
	stickyTumorInterfaceForces();
	shapeForces2D();
	adipocyteECMAdhesionForces();
}




/******************************

	T U M O R

	F I R E

*******************************/


void tumor2D::tumorFIRE(tumor2DMemFn forceCall, double Ftol, double dt0) {
	// local variables
	int i;

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
			else if (x[i] < 0.0 && pbc[i % NDIM])
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
		cout << endl << endl;
	}
}



/******************************

	P R O T O C O L S

*******************************/


void tumor2D::setupCheck(){
	// check NVTOT set up
	if (NVTOT <= 0){
		cerr << "** ERROR: in setupCheck, NVTOT = " << NVTOT << ". Ending. " << endl;
		exit(1);
	}

	// check tN
	if (tN > NCELLS){
		cerr << "** ERROR: in setupCheck, tN = " << tN << ", which is > NCELLS = " << NCELLS << ". Ending. " << endl;
		exit(1);
	}


	// check initialization of shape parameters
	if (l0.size() != NVTOT){
		cerr << "** ERROR: in setupCheck, l0.size = " << l0.size() << ", which is != NVTOT = " << NVTOT << ". Ending. " << endl;
		exit(1);
	}
	if (l0_init.size() != NVTOT){
		cerr << "** ERROR: in setupCheck, l0_init.size = " << l0_init.size() << ", which is != NVTOT = " << NVTOT << ". Ending. " << endl;
		exit(1);
	}
	if (t0.size() != NVTOT){
		cerr << "** ERROR: in setupCheck, t0.size = " << t0.size() << ", which is != NVTOT = " << NVTOT << ". Ending. " << endl;
		exit(1);
	}
	if (a0.size() != NCELLS){
		cerr << "** ERROR: in setupCheck, a0.size = " << a0.size() << ", which is != NCELLS = " << NCELLS << ". Ending. " << endl;
		exit(1);
	}
}


// compression with boundary forces
void tumor2D::tumorCompression(double Ftol, double Ptol, double dt0, double dphi0){
	// check correct setup
	setupCheck();

	// local variables
	int k = 0, nr;
	double pcheck=0.0, phi0, scaleFactor = 1.0;

	// initialize preferred packing fraction
	phi0 = vertexPreferredPackingFraction2D();

	// loop until pcheck > Ptol is found
	while (pcheck < Ptol && k < itmax) {
		// relax configuration (pass repsulive force update member function)
		tumorFIRE(&tumor2D::repulsiveTumorInterfaceForceUpdate, Ftol, dt0);

		// update pressure
		pcheck = 0.5 * (stress[0] + stress[1]);

		// remove rattlers
		nr = removeRattlers();

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << " 	T U M O R    2 D  							" << endl;
		cout << " 	  	I S O T R O P I C 						" << endl;
		cout << "			C O M P R E S S I O N 				" << endl << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi0 			= " << phi0 << endl;
		cout << "	* phi 			= " << vertexPackingFraction2D() << endl;
		cout << "	* scaleFactor 	= " << scaleFactor << endl;
		cout << "	* pcheck 		= " << pcheck << endl;
		cout << "	* U 		 	= " << U << endl;
		cout << "	* Nvv  			= " << vvContacts() << endl;
		cout << "	* Ncc 			= " << ccContacts() << endl;
		cout << endl << endl;

		// scale particle sizes
		scaleFactor = sqrt((phi0 + dphi0) / phi0);
		scaleParticleSizes2D(scaleFactor);

		// update packing fraction
		phi0 = vertexPreferredPackingFraction2D();

		// update iterate
		k++;
	}
}


// invasion protocol
// DEBUG FORCE CALL FROM interfaceInvasion.cpp
void tumor2D::invasion(tumor2DMemFn forceCall, double dDr, double dPsi, double Drmin, int NT, int NPRINTSKIP){
	// check correct setup
	setupCheck();

	// local variables
	int k, i, ci, cj;
	double t = 0.0, zta, Drtmp;

	// attach pins
	updateECMAttachments(1);

	// loop over time, have active brownian crawlers invade adipocytes
	for (k=0; k<NT; k++){
		// pbcs and reset forces
		for (i=0; i<vertDOF; i++){
			// recenter in box (only if y)
			if (i % NDIM == 1){
				if (x[i] > L[1])
					x[i] -= L[1];
				else if (x[i] < 0)
					x[i] += L[1];
			}
		}

		// update forces
		CALL_MEMBER_FN(*this, forceCall)();

		// update active brownian crawler
		activeBrownianCrawlerUpdate();

		// update positions (EULER UPDATE, OVERDAMPED)
		for (i=0; i<vertDOF; i++)
			x[i] += dt * F[i];

		// increase persistence + drift director if close to adipocytes
		for (ci=0; ci<tN; ci++){
			// get number of tumor-adipocyte contacts
			zta = 0.0;
			for (cj=tN; cj<NCELLS; cj++)
				zta += cij[NCELLS * ci + cj - (ci + 1) * (ci + 2) / 2]++;

			// change persistence
			Drtmp = Dr0*(1 - (zta/nv[ci])*dDr);
			if (Drtmp > Drmin)
				Dr[ci] = Drtmp;
			else
				Dr[ci] = Drmin;

			// change psi
			psi[ci] -= dt * (zta/nv[ci]) * dPsi * psi[ci];
		}

		// update psi based on persistence
		psiDiffusion();

		// update time
		t += dt;

		// print message console, print position to file
		if (k % NPRINTSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "			invading tumor cells 			" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** k 			= " << k << endl;
			cout << "	** p 			= " << 0.5*(stress[0] + stress[1]) << endl;
			cout << "	** phi 			= " << vertexPackingFraction2D() << endl;

			// print vertex positions to check placement
			cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
			printTumorInterface(t);
		}
	}
}


// just get them bois crawlin'
void tumor2D::crawling(tumor2DMemFn forceCall, int NT, int NPRINTSKIP){
	// check correct setup
	setupCheck();

	// local variables
	int k, i, ci, cj;
	double t = 0.0;

	// loop over time, have active brownian crawlers invade adipocytes
	for (k=0; k<NT; k++){
		// pbcs
		for (i=0; i<vertDOF; i++){
			if (x[i] > L[i % NDIM] && pbc[i % NDIM])
				x[i] -= L[i % NDIM];
			else if (x[i] < 0.0 && pbc[i % NDIM])
				x[i] += L[i % NDIM];
		}

		// update forces
		CALL_MEMBER_FN(*this, forceCall)();

		// update active brownian crawler
		activeBrownianCrawlerUpdate();

		// update positions (EULER UPDATE, OVERDAMPED)
		for (i=0; i<vertDOF; i++)
			x[i] += dt * F[i];

		// update psi based on persistence
		psiDiffusion();

		// update time
		t += dt;

		// print message console, print position to file
		if (k % NPRINTSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "			active tumor cells 				" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** k 			= " << k << endl;
			cout << "	** p 			= " << 0.5*(stress[0] + stress[1]) << endl;
			cout << "	** phi 			= " << vertexPackingFraction2D() << endl;

			// print vertex positions to check placement
			cout << "\t** PRINTING POSITIONS TO FILE... " << endl;
			printTumorCells(t);
		}
	}
}





/******************************

	P R I N T I N G

	F U N C T I O N S

*******************************/


void tumor2D::printTumorInterface(double t){
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
	posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << setw(w) << left << tN << endl;
	posout << setw(w) << left << "TSTEP" << setw(wnum) << setprecision(pnum) << left << dt << endl;
	posout << setw(w) << left << "TCURR" << setw(wnum) << setprecision(pnum) << left << t << endl;
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

	// print wall stress info
	posout << setw(w) << left << "WPRSS";
	posout << setw(wnum) << setprecision(pnum) << left << wpress.at(0);
	posout << setw(wnum) << setprecision(pnum) << left << wpress.at(1);
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
		if (ci < tN){
			posout << setw(wnum) << left << psi.at(ci);
			posout << setw(wnum) << left << Dr.at(ci);
		}
		else{
			posout << setw(wnum) << left << pinpos[NDIM*(ci-tN)];
			posout << setw(wnum) << left << pinpos[NDIM*(ci-tN) + 1];
			posout << setw(w) << left << pinattach[ci-tN];
		}
		posout << endl;

		// get initial vertex positions
		gi = gindex(ci, 0);
		xi = x.at(NDIM * gi);
		yi = x.at(NDIM * gi + 1);

		// place back in box center
		if (pbc[0])
			xi = fmod(xi, Lx);
		if (pbc[1])
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
	posout << setw(w) << left << "ENDFR" << " " << endl;
}

void tumor2D::printTumorCells(double t){
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
	posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << setw(w) << left << tN << endl;
	posout << setw(w) << left << "TSTEP" << setw(wnum) << setprecision(pnum) << left << dt << endl;
	posout << setw(w) << left << "TCURR" << setw(wnum) << setprecision(pnum) << left << t << endl;
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
		posout << setw(wnum) << left << psi.at(ci);
		posout << setw(wnum) << left << Dr.at(ci);
		posout << endl;

		// get initial vertex positions
		gi = gindex(ci, 0);
		xi = x.at(NDIM * gi);
		yi = x.at(NDIM * gi + 1);

		// place back in box center
		if (pbc[0])
			xi = fmod(xi, Lx);
		if (pbc[1])
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
	posout << setw(w) << left << "ENDFR" << " " << endl;
}


