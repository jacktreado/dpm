/*

	FUNCTION DEFINITIONS for MESO class

	Jack Treado, 06/09/21

*/

#include "meso2D.h"

// namespace
using namespace Eigen;
using namespace std;


/******************************

	R E A D - I N

	C O N S T R U C T O R

*******************************/

meso2D::meso2D(string &inputFileStr,int seed) : dpm(2) {
	// open file
	ifstream inputobj(inputFileStr.c_str());
	if (!inputobj.is_open()){
		cerr << "** ERROR: In meso2D constructor, could not open file " << inputFileStr << ", ending here. " << endl;
		exit(1);
	}

	// set variables to default
	NBUBBLES = 0;
	betaEff=0.0; ctcdel=1.0; ctch=0.5; cL=0.0; aL=1.0; cB=0.0; cKb=0.0; pbc[0]=1; pbc[1]=1;
	ka = 1.0;
	kl = 1.0;
	kb = 0.0;
	kc = 1.0;

	// local variables
	int nvtmp, ci, vi, i;
	double val;
	double lxtmp, lytmp;
	double s1, s2, s3;
	double wp1, wp2;
	double a0tmp, l0tmp, t0tmp, kbitmp;
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

	// verify input file
	if (NCELLS < 1){
		cerr << "** ERROR: in meso2D constructor, NCELLStmp = " << NCELLS << ". Ending here." << endl;
		exit(1);
	}
	
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
	Pinst = 0.0;

	// szList and nv (keep track of global vertex indices)
	nv.resize(NCELLS);
	szList.resize(NCELLS);
	a0.resize(NCELLS);

	// initialize NVTOT to 0
	NVTOT = 0;

	// loop over cells, read in coordinates
	cout << "\n\n** LOOPING OVER DATA, PRINTING INFO..." << endl;
	for (ci=0; ci<NCELLS; ci++){
		// first parse cell info
		getline(inputobj, inputStr);
		sscanf(inputStr.c_str(),"CINFO %d %*d %lf %*lf %*lf",&nvtmp,&a0tmp);

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
			sscanf(inputStr.c_str(),"VINFO %*d %*d %lf %lf %lf %lf %lf %lf %*lf",&xtmp,&ytmp,&rtmp,&l0tmp,&t0tmp,&kbitmp);

			// print to console
			cout << "\t ** " << inputStr << endl;

			// push back 
			x.push_back(xtmp);
			x.push_back(ytmp);
			r.push_back(rtmp);
			l0.push_back(l0tmp);
			t0.push_back(t0tmp);
			kbi.push_back(kbitmp);
		}
	}
	vertDOF = NDIM * NVTOT;
	cout << "** NVTOT = " << NVTOT << ", vertDOF = " << vertDOF << ", NCELLS = " << NCELLS << endl;
	cout << "** FINISHED LOOPING OVER DATA, CLOSING INPUT FILE OBJ\n" << endl;

	// initialize vertex indexing
	initializeVertexIndexing2D();

	// initialize contact information
	cij.resize(NCELLS * (NCELLS - 1) / 2);
	zv.resize(NVTOT);
	zc.resize(NVTOT);
	
	fill(cij.begin(),cij.end(),0);
	fill(zv.begin(), zv.end(), 0);
	fill(zc.begin(), zc.end(), 0);

	// initialize contact network for bonds
	initializeVertexContactNetwork();

	// close input file object
	inputobj.close();

	// seed random number generator
	srand48(seed);
}

meso2D::meso2D(string &inputFileStr, double voidBubbleScale, int seed) : dpm(2) {
	// open file
	ifstream inputobj(inputFileStr.c_str());
	if (!inputobj.is_open()){
		cerr << "** ERROR: In meso2D void Bubble constructor, could not open file " << inputFileStr << ", ending here. " << endl;
		exit(1);
	}

	// set variables to default
	NBUBBLES = 0;
	betaEff=0.0; ctcdel=1.0; ctch=0.5; cL=0.0; aL=1.0; cB=0.0; cKb=0.0; pbc[0]=1; pbc[1]=1;
	ka = 1.0;
	kl = 1.0;
	kb = 0.0;
	kc = 1.0;

	// local variables
	int nctmp, nvtmp, ci, vi, i;
	double val;
	double lxtmp, lytmp;
	double s1, s2, s3;
	double wp1, wp2;
	double a0tmp, l0tmp, t0tmp, kbitmp;
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
	sscanf(inputStr.c_str(),"NUMCL %d",&nctmp);
	cout << "\t ** " << inputStr << endl;

	getline(inputobj, inputStr);
	sscanf(inputStr.c_str(),"NUMBL %d",&NBUBBLES);
	cout << "\t ** " << inputStr << endl;

	// let number of cells be # true cells + # bubbles
	NCELLS = nctmp + NBUBBLES;

	// verify input file
	if (NCELLS < 1){
		cerr << "** ERROR: in meso2D constructor, NCELLStmp = " << NCELLS << ". Ending here." << endl;
		exit(1);
	}
	
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
	Pinst = 0.0;

	// szList and nv (keep track of global vertex indices)
	nv.resize(NCELLS);
	szList.resize(NCELLS);
	a0.resize(NCELLS);

	// initialize NVTOT to 0
	NVTOT = 0;

	// loop over cells, read in coordinates
	cout << "\n\n** LOOPING OVER DATA, PRINTING INFO..." << endl;
	for (ci=0; ci<NCELLS; ci++){
		// first parse cell info
		getline(inputobj, inputStr);
		sscanf(inputStr.c_str(),"CINFO %d %*d %lf %*lf %*lf",&nvtmp,&a0tmp);

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
			sscanf(inputStr.c_str(),"VINFO %*d %*d %lf %lf %lf %lf %lf %lf %*lf",&xtmp,&ytmp,&rtmp,&l0tmp,&t0tmp,&kbitmp);

			// print to console
			cout << "\t ** " << inputStr << endl;

			// push back 
			x.push_back(xtmp);
			x.push_back(ytmp);
			r.push_back(rtmp);
			l0.push_back(l0tmp);
			t0.push_back(t0tmp);
			kbi.push_back(kbitmp);
		}
	}
	vertDOF = NDIM * NVTOT;
	cout << "** NVTOT = " << NVTOT << ", vertDOF = " << vertDOF << ", NCELLS = " << NCELLS << ", NBUBBLES = " << NBUBBLES << ", NTRUECELLS = " << nctmp << endl;
	cout << "** FINISHED LOOPING OVER DATA, CLOSING INPUT FILE OBJ\n" << endl;

	// initialize vertex indexing
	initializeVertexIndexing2D();

	// initialize contact information
	cij.resize(NCELLS * (NCELLS - 1) / 2);
	zv.resize(NVTOT);
	zc.resize(NVTOT);
	
	fill(cij.begin(),cij.end(),0);
	fill(zv.begin(), zv.end(), 0);
	fill(zc.begin(), zc.end(), 0);

	// initialize contact network for bonds
	initializeVertexContactNetwork();

	// close input file object
	inputobj.close();

	// seed random number generator
	srand48(seed);
}



// overloaded equal operator (copy everything except outputs) (for Lees-Edwards shear test)
void meso2D::operator=(const meso2D &rhs){
	// -- DPM 2D VARIABLES

	// dpm int scalars
	NCELLS = rhs.NCELLS;
	NDIM = rhs.NDIM;
	NNN = rhs.NNN;
	NVTOT = rhs.NVTOT;
	vertDOF = rhs.vertDOF;

	// time step size
	dt = rhs.dt;

	// copy energetic info
	U = rhs.U;
	ka = rhs.ka;
	kl = rhs.kl;
	kb = rhs.kb;
	kc = rhs.kc;
	l1 = rhs.l1;
	l2 = rhs.l2;

	// boundary parameters
	L = rhs.L;
	pbc = rhs.pbc;

	// particle shape parameters
	a0 = rhs.a0;
	l0 = rhs.l0;
	t0 = rhs.t0;
	r = rhs.r;

	// indexing variables
	nv = rhs.nv;
	szList = rhs.szList;
	im1 = rhs.im1;
	ip1 = rhs.ip1;

	// dynamical variables
	x = rhs.x;
	v = rhs.v;
	F = rhs.F;

	// stress
	stress = rhs.stress;

	// cell-cell contact network
	cij = rhs.cij;

	// box linked-list variables
	NBX = rhs.NBX;
	sb = rhs.sb;
	lb = rhs.lb;
	nn = rhs.nn;
	head = rhs.head;
	last = rhs.last;
	list = rhs.list;

	// -- MESO 2D VARIABLES

	// max number of vertices
	NVMAX = rhs.NVMAX;

	// bending energy per vertex
	kbi = rhs.kbi;

	// vertex-vertex bond network
	gij = rhs.gij;

	// vv contacts per cell
	zc = rhs.zc;
	zv = rhs.zv;

	// adhesion parameters
	betaEff = rhs.betaEff;
	ctcdel = rhs.ctcdel;
	ctch = rhs.ctch;

	// aging during development
	cL = rhs.cL;
	aL = rhs.aL;
	cB = rhs.cB;
	cKb = rhs.cKb;
}





/******************************

	I N I T I A L -

	I Z A T I O N

*******************************/

void meso2D::initializeMesophyllCells(double dispersion, double calA0, double phi0, double Ftol, int n1){
	// print message to console
	cout << "	-- IN MESO2D FUNCTION initializeMesophyllCells, dispersion = " << dispersion << ", calA0 = " << calA0 << ", phi0 = " << phi0 << ", Ftol = " << Ftol << ", n1 = " << n1 << endl;

	// set initial mechanical constants
	ka = 1.0;
	kl = 1.0;
	kb = 0.0;
	kc = 1.0;

	// initialize particles drawn from Gaussian distribution
	gaussian2D(dispersion, calA0, n1);

	// initialize contact network vector
	initializeVertexContactNetwork();

	// initialize cell centers using SP FIRE
	initializePositions2D(phi0, Ftol);

	// set list of bending energies to 0
	kbi.resize(NVTOT);
	setkbi(0.0);

	// initialize vertex-vertex contact list
	zv.resize(NVTOT);
	fill(zv.begin(), zv.end(), 0);

	// just in case, initialize pressure to 0
	Pinst = 0.0;
}


void meso2D::initializeVertexContactNetwork(){
	// local variables
	int NVCTS, gi;

	// check that vertDOF has been assigned
	if (NVTOT <= 0){
		cerr << "	** ERROR: in meso2D::initializeVertexContactNetwork, NVTOT not assigned. Ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cerr << "	** ERROR: in meso2D::initializeVertexContactNetwork, vertDOF not assigned. Ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cerr << "	** ERROR: in meso2D::initializeVertexContactNetwork, nv vector not assigned. Ending here." << endl;
		exit(1);
	}

	// number of vertex-vertex contacts
	NVCTS = NVTOT*(NVTOT-1)/2;

	// initialize contact network
	gij.resize(NVCTS);
	for (gi=0; gi<NVCTS; gi++)
		gij.at(gi) = 0;
}







/******************************

	M E S O P H Y L L

	E D I T O R S  &

	U P D A T E S

*******************************/

double meso2D::meanl0(){
	int gi;
	double val;

	val = 0.0;
	for (gi=0; gi<NVTOT; gi++)
		val += l0[gi];
	val /= NVTOT;

	return val;
}

double meso2D::meancalA0(){
	int gi, ci, vi;
	double calA0tmp, val;

	val = 0.0;
	for (ci=0; ci<NCELLS; ci++){
		calA0tmp = 0.0;
		for (vi=0; vi<nv[ci]; vi++)
			calA0tmp += l0[szList[ci]+vi];
		calA0tmp *= calA0tmp;
		calA0tmp /= 4.0*PI*a0[ci];
		val += calA0tmp;
	}
	val /= NCELLS;

	return val;
}

double meso2D::meant0(){
	int gi;
	double val;

	val = 0.0;
	for (gi=0; gi<NVTOT; gi++)
		val += t0[gi];
	val /= NVTOT;

	return val;
}

double meso2D::meankb(){
	int gi;
	double val;

	val = 0.0;
	for (gi=0; gi<NVTOT; gi++)
		val += kbi[gi];
	val /= NVTOT;

	return val;
}








/******************************

	M E S O P H Y L L

	C E L L

	I N T E R A C T I O N S

*******************************/


// initialize bonds between overlapping vertices
void meso2D::initializeMesophyllBondNetwork(){
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, dx, dy;

	// sort particles
	sortNeighborLinkedList2D();

	// reset vertex-vertex bond network
	fill(gij.begin(), gij.end(), 0);
	fill(zc.begin(), zc.end(), 0);
	fill(zv.begin(), zv.end(), 0);

	// loop over boxes in cell linked list
	for (bi=0; bi<NBX; bi++){

		// get start of list of vertices
		pi = head[bi];

		// loop over linked list
		while (pi > 0){
			// real particle index
			gi = pi - 1;

			// next particle in list
			pj = list[pi];

			// loop down neighbors of pi in same interaction cell
			while (pj > 0){
				// real index of pj
				gj = pj - 1;

				// get cell and vert indices
				cindices(ci,vi,gi);
				cindices(cj,vj,gj);

				if (ci == cj){
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

							// add bond to network
							if (gi > gj)
								gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
							else if (gi < gj)
								gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1; 

							// update cell and vertex contacts
							zc[ci]++;
							zc[cj]++;

							zv[gi]++;
							zv[gj]++;
						}
					}
				}
				// update pj
				pj = list[pj];
			}

			// test overlaps with forward neighboring interaction cells
			for (bj=0; bj<NNN; bj++){
				// only check if boundaries permit
				if (nn[bi][bj] == -1)
					continue;

				// get first particle in neighboring interaction cell
				pj = head[nn[bi][bj]];

				// loop down neighbors of pi in same interaction cell
				while (pj > 0){
					// real index of pj
					gj = pj - 1;

					// get cell and vert indices
					cindices(ci,vi,gi);
					cindices(cj,vj,gj);

					if (ci == cj){
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
								// add bond
								if (gi > gj)
									gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
								else if (gi < gj)
									gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1; 

								// update cell and vertex contacts
								zc[ci]++;
								zc[cj]++;

								zv[gi]++;
								zv[gj]++;
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
}


// initialize bonds between overlapping CELLS
void meso2D::initializeMesoBubbleBondNetwork(){
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, dx, dy;
	int NTCELLS = NCELLS - NBUBBLES;

	// sort particles
	sortNeighborLinkedList2D();

	// reset vertex-vertex bond network
	fill(gij.begin(), gij.end(), 0);
	fill(zc.begin(), zc.end(), 0);
	fill(zv.begin(), zv.end(), 0);

	// loop over boxes in cell linked list
	for (bi=0; bi<NBX; bi++){

		// get start of list of vertices
		pi = head[bi];

		// loop over linked list
		while (pi > 0){
			// real particle index
			gi = pi - 1;

			// next particle in list
			pj = list[pi];

			// loop down neighbors of pi in same interaction cell
			while (pj > 0){
				// real index of pj
				gj = pj - 1;

				// get cell and vert indices
				cindices(ci,vi,gi);
				cindices(cj,vj,gj);

				if (ci == cj || ci >= NTCELLS || cj >= NTCELLS){
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

							// add bond to network
							if (gi > gj)
								gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
							else if (gi < gj)
								gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1; 

							// update cell and vertex contacts
							zc[ci]++;
							zc[cj]++;

							zv[gi]++;
							zv[gj]++;
						}
					}
				}
				// update pj
				pj = list[pj];
			}

			// test overlaps with forward neighboring interaction cells
			for (bj=0; bj<NNN; bj++){
				// only check if boundaries permit
				if (nn[bi][bj] == -1)
					continue;

				// get first particle in neighboring interaction cell
				pj = head[nn[bi][bj]];

				// loop down neighbors of pi in same interaction cell
				while (pj > 0){
					// real index of pj
					gj = pj - 1;

					// get cell and vert indices
					cindices(ci,vi,gi);
					cindices(cj,vj,gj);

					if (ci == cj || ci >= NTCELLS || cj >= NTCELLS){
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
								// add bond
								if (gi > gj)
									gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
								else if (gi < gj)
									gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1; 

								// update cell and vertex contacts
								zc[ci]++;
								zc[cj]++;

								zv[gi]++;
								zv[gj]++;
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
}

// mesophyll repulsive vertex forces
void meso2D::mesoRepulsiveVertexForces(){
	// local variables
	int ci, cj, gi, gj, vi, vj, bi, bj, pi, pj, boxid, sbtmp;
	double sij, rij, dx, dy;
	double ftmp, fx, fy;

	// sort particles
	sortNeighborLinkedList2D();

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
							ftmp = (kc / sij) * (1 - (rij / sij));
							fx = ftmp * (dx / rij);
							fy = ftmp * (dy / rij);

							// add to forces
							F[NDIM * gi] -= fx;
							F[NDIM * gi + 1] -= fy;

							F[NDIM * gj] += fx;
							F[NDIM * gj + 1] += fy;

							// increase potential energy
							U += 0.5 * kc * pow((1 - (rij / sij)), 2.0);

							// add to dUdL contribution to pressure
							Pinst -= ftmp*(rij/L[0]);

							// add to virial stress
							stress[0] += (dx * fx)/(L[0] * L[1]);
							stress[1] += (dy * fy)/(L[0] * L[1]);;
							stress[2] += (0.5 * (dx * fy + dy * fx))/(L[0] * L[1]);

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
								ftmp = (kc / sij) * (1 - (rij / sij));
								fx = ftmp * (dx / rij);
								fy = ftmp * (dy / rij);

								// add to forces
								F[NDIM * gi] -= fx;
								F[NDIM * gi + 1] -= fy;

								F[NDIM * gj] += fx;
								F[NDIM * gj + 1] += fy;

								// increae potential energy
								U += 0.5 * kc * pow((1 - (rij / sij)), 2.0);

								// add to dUdL contribution to pressure
								Pinst -= ftmp*(rij/L[0]);

								// add to virial stress
								stress[0] += (dx * fx)/(L[0] * L[1]);
								stress[1] += (dy * fy)/(L[0] * L[1]);;
								stress[2] += (0.5 * (dx * fy + dy * fx))/(L[0] * L[1]);

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
}

// mesophyll specific shape forces
void meso2D::mesoShapeForces(){
	// local variables
	int ci, gi, vi, nvtmp;
	double fa, fli, flim1, fbi, fbim1, cx, cy, xi, yi;
	double flx, fly, fbx, fby;
	double l0im1, l0i, a0tmp, atmp;
	double dx, dy, da, dli, dlim1, dtim1, dti, dtip1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
	double nim1x, nim1y, nix, niy, sinim1, sini, sinip1, cosim1, cosi, cosip1;
	double ddtim1, ddti;

	// loop over vertices, add to force
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
				U += 0.5 * ka * (da * da);

				// update pressure
				Pinst += ((2.0 * ka * atmp)/(a0tmp * L[0])) * da;

				// shape force parameters
				fa = ka * (da / a0tmp);

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
		flim1 	= kl/l0im1;
		fli 	= kl/l0i;

		// add to forces
		flx 			= (fli*dli*lix/li) - (flim1*dlim1*lim1x/lim1);
		fly 			= (fli*dli*liy/li) - (flim1*dlim1*lim1y/lim1);
		F[NDIM*gi] 		+= flx;
		F[NDIM*gi + 1] 	+= fly;
		
		// update potential energy
		U += 0.5 * kl *(dli * dli);

		// update pressure
		Pinst += ((kl * li)/(L[0] * l0i)) * dli;

		// -- Bending force
		fbi = kbi[gi];
		fbim1 = kbi[im1[gi]];

		if (fbi > 0 || fbim1 > 0){
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
			sinim1 = lim2x*lim1y - lim2y*lim1x;
			cosim1 = lim2x*lim1x + lim2y*lim1y;

			sini = lim1x*liy - lim1y*lix;
			cosi = lim1x*lix + lim1y*liy;

			sinip1 = lix*lip1y - liy*lip1x;
			cosip1 = lix*lip1x + liy*lip1y;

			// get normal vectors
			nim1x = lim1y/lim1;
			nim1y = -lim1x/lim1;

			nix = liy/li;
			niy = -lix/li;

			// get change in angles
			dtim1 = atan2(sinim1,cosim1) - t0[im1[gi]];
			dti = atan2(sini,cosi) - t0[gi];
			dtip1 = atan2(sinip1,cosip1) - t0[ip1[gi]];

			// get delta delta theta's
			ddtim1 = (dti - dtim1)/lim1;
			ddti = (dtip1 - dti)/li;

			// add to force
			fbx 			= fbi*ddti*nix - fbim1*ddtim1*nim1x;
			fby 			= fbi*ddti*niy - fbim1*ddtim1*nim1y;
			F[NDIM*gi] 		+= fbx;
			F[NDIM*gi + 1] 	+= fby;

			// update potential energy
			U += 0.5 * kbi[gi] * (dti * dti);
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

// mesophyll specific shape forces WITH APPLIED SHEAR STRAIN GAMMA (assume pbcs)
// NOTE: add Pinst contribution to pressure measurement, need to figure out Sinst (shear stress) from shape contribution
void meso2D::mesoShapeForces(double gamma){
	// local variables
	int ci, gi, vi, nvtmp, im;
	double fa, fli, flim1, fbi, fbim1, cx, cy, xi, yi;
	double flx, fly, fbx, fby;
	double l0im1, l0i, a0tmp, atmp;
	double dx, dy, da, dli, dlim1, dtim1, dti, dtip1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
	double nim1x, nim1y, nix, niy, sinim1, sini, sinip1, cosim1, cosi, cosip1;
	double ddtim1, ddti;

	// loop over vertices, add to force
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
				atmp = area(ci,gamma);
				da = (atmp/a0tmp) - 1.0;

				// update potential energy
				U += 0.5 * ka * (da * da);

				// shape force parameters
				fa = ka * (da / a0tmp);

				// compute cell center of mass
				xi = x[NDIM*gi];
				yi = x[NDIM*gi + 1];
				cx = xi; 
				cy = yi;
				for (vi=1; vi<nvtmp; vi++){
					// get distances between vim1 and vi
					dy = x[NDIM*(gi+vi) + 1] - yi;
					im = round(dy/L[1]);
					dy -= L[1]*im;

					dx = x[NDIM*(gi+vi)] - xi;
					dx -= L[1]*im*gamma;
					dx -= L[0]*round(dx/L[0]);

					// add to centers
					xi += dx;
					yi += dy;

					cx += xi;
					cy += yi;
				}
				cx /= nvtmp;
				cy /= nvtmp;

				// get coordinates relative to center of mass
				riy = x[NDIM*gi + 1] - cy;
				im = round(riy/L[1]);
				riy -= L[1]*im;

				rix = x[NDIM*gi] - cx;
				rix -= L[1]*im*gamma;
				rix -= L[0]*round(rix/L[0]);



				// get prior adjacent vertices
				rim2y = x[NDIM*im1[im1[gi]] + 1] - cy;
				im = round(rim2y/L[1]);
				rim2y -= L[1]*im;

				rim2x = x[NDIM*im1[im1[gi]]] - cx;
				rim2x -= L[1]*im*gamma;
				rim2x -= L[0]*round(rim2x/L[0]);


				rim1y = x[NDIM*im1[gi] + 1] - cy;
				im = round(rim1y/L[1]);
				rim1y -= L[1]*im;

				rim1x = x[NDIM*im1[gi]] - cx;
				rim1x -= L[1]*im*gamma;
				rim1x -= L[0]*round(rim1x/L[0]);



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
		rip1y = x[NDIM*ip1[gi] + 1] - cy;
		im = round(rip1y/L[1]);
		rip1y -= L[1]*im;

		rip1x = x[NDIM*ip1[gi]] - cx;
		rip1x -= L[1]*im*gamma;
		rip1x -= L[0]*round(rip1x/L[0]);

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
		flim1 	= kl/l0im1;
		fli 	= kl/l0i;

		// add to forces
		flx 			= (fli*dli*lix/li) - (flim1*dlim1*lim1x/lim1);
		fly 			= (fli*dli*liy/li) - (flim1*dlim1*lim1y/lim1);
		F[NDIM*gi] 		+= flx;
		F[NDIM*gi + 1] 	+= fly;

		// update potential energy
		U += 0.5 * kl * (dli * dli);

		// -- Bending force
		fbi = kbi[gi];
		fbim1 = kbi[im1[gi]];

		if (fbi > 0 || fbim1 > 0){
			// get ip2 for third angle
			rip2y = x[NDIM*ip1[ip1[gi]] + 1] - cy;
			im = round(rip2y/L[1]);
			rip2y -= L[1]*im;

			rip2x = x[NDIM*ip1[ip1[gi]]] - cx;
			rip2x -= L[1]*im*gamma;
			rip2x -= L[0]*round(rip2x/L[0]);


			// get last segment length
			lip1x = rip2x - rip1x;
			lip1y = rip2y - rip1y;

			// get angles
			sinim1 = lim2x*lim1y - lim2y*lim1x;
			cosim1 = lim2x*lim1x + lim2y*lim1y;

			sini = lim1x*liy - lim1y*lix;
			cosi = lim1x*lix + lim1y*liy;

			sinip1 = lix*lip1y - liy*lip1x;
			cosip1 = lix*lip1x + liy*lip1y;

			// get normal vectors
			nim1x = lim1y/lim1;
			nim1y = -lim1x/lim1;

			nix = liy/li;
			niy = -lix/li;

			// get change in angles
			dtim1 = atan2(sinim1,cosim1) - t0[im1[gi]];
			dti = atan2(sini,cosi) - t0[gi];
			dtip1 = atan2(sinip1,cosip1) - t0[ip1[gi]];

			// get delta delta theta's
			ddtim1 = (dti - dtim1)/lim1;
			ddti = (dtip1 - dti)/li;

			// add to force
			fbx 			= fbi*ddti*nix - fbim1*ddtim1*nim1x;
			fby 			= fbi*ddti*niy - fbim1*ddtim1*nim1y;
			F[NDIM*gi] 		+= fbx;
			F[NDIM*gi + 1] 	+= fby;

			// update potential energy
			U += 0.5 * kbi[gi] * (dti * dti);
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


// update forces for bonded mesophyll cells
// * Shape forces of individual cells
// * Repulsive interactions between overlapping vertices
// * Interactions between bonded vertices
void meso2D::mesoNetworkForceUpdate(){
	// local variables
	int i, gi, gj, ci, cj, vi, vj;
	double rij, sij, zij, dx, dy, fx, fy, ftmp, Fvirial;

	// normal update (shape + repulsive forces) from base class
	resetForcesAndEnergy();
	Pinst = 0.0;
	mesoShapeForces();
	mesoRepulsiveVertexForces();

	// update bonded forces
	for (gi=0; gi<NVTOT; gi++){
		for (gj=gi+1; gj<NVTOT; gj++){
			if (gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2]){

				// contact distance
				sij = r[gi] + r[gj];

				// get vertex-vertex distance
				dx = x[NDIM*gj] - x[NDIM*gi];
				if (pbc[0])
					dx -= L[0]*round(dx/L[0]);

				dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
				if (pbc[1])
					dy -= L[1]*round(dy/L[1]);

				rij = sqrt(dx*dx + dy*dy);

				// only compute force if spring is extended
				if (rij > sij){
					// get cell indices
					cindices(ci,vi,gi);
					cindices(cj,vj,gj);

					// zij: determines strength of bond attraction
					zij = 0.5*(zc[ci] + zc[cj])*ctcdel + 1.0;

					// force scale
					ftmp 				= (kc/(sij * zij))*(1 - (rij/sij));
					fx 					= ftmp*(dx/rij);
					fy 					= ftmp*(dy/rij);

					// add to forces
					F[NDIM*gi] 			-= fx;
					F[NDIM*gi + 1] 		-= fy;

					F[NDIM*gj] 			+= fx;
					F[NDIM*gj + 1] 		+= fy;

					// add to cell-cell contacts
					if (ci > cj)
						cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]++;
					else if (ci < cj)
						cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]++; 

					// increae potential energy
					U += 0.5*(kc/zij)*pow((1 - (rij/sij)),2.0);

					// add to virial stress tensor
					stress[0] += (dx*fx)/(L[0]*L[1]);
					stress[1] += (dy*fy)/(L[0]*L[1]);
					stress[2] += (0.5*(dx*fy + dy*fx))/(L[0]*L[1]);

					// add to dUdL
					Pinst -= ftmp*(rij/L[0]);
				}
			}
		}
	}

	// virial force contribution to pressure
	Fvirial = 0.0;
	for (i=0; i<vertDOF; i++)
		Fvirial += F[i]*x[i];
	Fvirial /= 2.0*L[0]*L[1];

	// finalize computation of pressure
	Pinst = Fvirial - (Pinst/(2.0*L[0]));
}


// update forces AT FIXED SHEAR STRAIN (assume Lees-Edwards boundary conditions, fixed contact network)
// NOTE: add Pinst contribution to pressure measurement, need to figure out Sinst (shear stress) from shape contribution
void meso2D::mesoNetworkForceUpdate(double gamma, vector<bool> &gijtmp){
	// local variables
	int gi, gj, ci, cj, vi, vj, im, NVVCTS;
	double rij, sij, zij, dx, dy, fx, fy, ftmp;

	// normal update (shape + repulsive forces) from base class
	resetForcesAndEnergy();
	mesoShapeForces(gamma);

	// update bonded forces
	for (gi=0; gi<NVTOT; gi++){
		for (gj=gi+1; gj<NVTOT; gj++){
			if (gijtmp[NVTOT*gi + gj - (gi+1)*(gi+2)/2] == 1){

				// get cell indices
				cindices(ci,vi,gi);
				cindices(cj,vj,gj);

				// contact distance
				sij = r[gi] + r[gj];

				// get vertex-vertex distance
				dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
				im = round(dy/L[1]);
				dy -= L[1]*im;

				dx = x[NDIM*gj] - x[NDIM*gi];
				dx -= L[1]*im*gamma;
				dx -= L[0]*round(dx/L[0]);

				// get true distance
				rij = sqrt(dx*dx + dy*dy);

				// compute forces and potential
				if (rij > sij){
					// zij: determines strength of bond attraction
					zij 	= 0.5*(zc[ci] + zc[cj])*ctcdel + 1.0;

					// force scale
					ftmp 	= (kc/(sij * zij))*(1 - (rij/sij));
					fx 		= ftmp*(dx/rij);
					fy 		= ftmp*(dy/rij);

					// increae potential energy
					U 		+= 0.5*(kc/zij)*pow((1 - (rij/sij)),2.0);
				}
				else {
					// forces
					ftmp 	= (kc/sij)*(1 - (rij/sij));
					fx 	 	= ftmp*(dx/rij);
					fy 		= ftmp*(dy/rij);

					// increae potential energy
					U 		+= 0.5*kc*pow((1 - (rij/sij)),2.0);
				}

				// add to total forces
				F[NDIM*gi] 			-= fx;
				F[NDIM*gi + 1] 		-= fy;

				F[NDIM*gj] 			+= fx;
				F[NDIM*gj + 1] 		+= fy;

				// add to virial stress tensor
				stress[0] += (dx*fx)/(L[0]*L[1]);
				stress[1] += (dy*fy)/(L[0]*L[1]);
				stress[2] += (0.5*(dx*fy + dy*fx))/(L[0]*L[1]);
			}
		}
	}
}


// update forces for mesophyll cells with pin force (ASSUME PIN ALREADY ASSIGNED)
void meso2D::mesoPinForceUpdate(vector<double>& xpin, double kcspring){
	// local variables
	int gi, ci, vi, nvtmp;
	double cx, cy, dcx, dcy, fx, fy;

	// update network forces
	mesoNetworkForceUpdate();

	// update sinking forces
	for (ci=0; ci<NCELLS; ci++){
		// tmp nv
		nvtmp = nv[ci];

		// get center of mass coordinates
		com2D(ci,cx,cy);

		// get distances to pins
		dcx = cx - xpin[NDIM*ci];
		if (pbc[0])
			dcx -= L[0]*round(dcx/L[0]);

		dcy = cy - xpin[NDIM*ci + 1];
		if (pbc[1])
			dcy -= L[1]*round(dcy/L[1]);

		// get forces
		fx = -kcspring*dcx;
		fy = -kcspring*dcy;

		// add to virial stress tensor
		stress[0] += (dcx*fx)/(L[0]*L[1]);
		stress[1] += (dcy*fy)/(L[0]*L[1]);
		stress[2] += (0.5*(dcx*fy + dcy*fx))/(L[0]*L[1]);


		// add to forces 
		gi = szList[ci];
		for (vi=0; vi<nvtmp; vi++){
			F[NDIM*gi] += fx/nvtmp;
			F[NDIM*gi + 1] += fy/nvtmp;
			gi++;
		}
	}
}



/******************************

	M E S O P H Y L L

	I N T E G R A T O R S

*******************************/


// use FIRE to minimize potential energy in mesophyll cell network
void meso2D::mesoFIRE(meso2DMemFn forceCall, double Ftol, double dt0){
	// local variables
	int i;

	// check to see if cell linked-list has been initialized
	if (NBX == -1){
		cerr << "	** ERROR: In meso2D::mesoNetworkFIRE, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
		exit(1);
	}

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

	dtmax   	= 20.0*dt;
	dtmin   	= 1e-2*dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;

	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);

	// relax forces using FIRE
	while ((fcheck > Ftol || fireit < NDELAY) && fireit < itmax){
		// compute P
		P = 0.0;
		for (i=0; i<vertDOF; i++)
			P += v[i]*F[i];

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "		M E S O P H Y L L 			" << endl;
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
		for (i=0; i<vertDOF; i++)
			x[i] += dt*v[i];

		// update forces (function passed as argument)
		CALL_MEMBER_FN(*this, forceCall)();

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
		cout << "		M E S O P H Y L L 			" << endl;
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

// use FIRE to minimize ENTHALPY, by varying box size as well as coordinates + void perimeters
void meso2D::mesoEnthalpyFIRE(meso2DMemFn forceCall, double Ftol, double dPtol, double P0, double dt0){
	// local variables
	int i, d;

	// box size, momentum, and internal virial pressure
	double V=L[0]*L[1], Pi=0.0, P=0.0;

	// check to see if cell linked-list has been initialized
	if (NBX == -1){
		cerr << "	** ERROR: In meso2D::mesoEnthalpyFIRE, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
		exit(1);
	}

	// FIRE variables
	double PFIRE, fnorm, fcheck, fcheckFrc, dPcheck, vnorm, alpha, dtmax, dtmin;
	int npPos, npNeg, fireit;

	// set dt based on geometric parameters
	setdt(dt0);

	// Initialize FIRE variables
	PFIRE  		= 0;	
	fnorm 		= 0;
	vnorm 		= 0;
	alpha   	= alpha0;

	dtmax   	= 20.0*dt;
	dtmin   	= 1e-1*dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;
	fcheckFrc 	= fcheck;

	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);
	Pinst = 0.0;
	CALL_MEMBER_FN(*this, forceCall)();
	// P = Pinst;
	P = 0.5*(stress[0] + stress[1]);

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
			cout << "		M E S O P H Y L L 					" << endl;
			cout << " 	F I R E 								" << endl;
			cout << " 		E N T H A L P Y  					" << endl;
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

		// update forces (function passed as argument)
		Pinst = 0.0;
		CALL_MEMBER_FN(*this, forceCall)();

		// update instantaneous pressure
		// P = 0.5*(stress[0] + stress[1]);
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
		cout << "		M E S O P H Y L L 					" << endl;
		cout << " 	F I R E 								" << endl;
		cout << " 		E N T H A L P Y  					" << endl;
		cout << "	M I N I M I Z A T I O N 				" << endl;
		cout << "		C O N V E R G E D! 					" << endl;
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


// FIRE minimization at fixed shear strain
void meso2D::mesoShearStrainEnthalpyFIRE(double gamma, double Ftol, double P0, double dt0, vector<bool> &gijtmp){
	// local variables
	int i, d;
	const double dPtol = 1e-10;

	// box size, momentum, and internal virial pressure
	double V=L[0]*L[1], Pi=0.0, P=0.0;

	// check to see if cell linked-list has been initialized
	if (NBX == -1){
		cerr << "	** ERROR: In meso2D::mesoShearStrainEnthalpyFIRE, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
		exit(1);
	}

	// FIRE variables
	double PFIRE, fnorm, fcheck, fcheckFrc, dPcheck, vnorm, alpha, dtmax, dtmin;
	int npPos, npNeg, fireit;

	// set dt based on geometric parameters
	setdt(dt0);

	// Initialize FIRE variables
	PFIRE  		= 0;	
	fnorm 		= 0;
	vnorm 		= 0;
	alpha   	= alpha0;

	dtmax   	= 20.0*dt;
	dtmin   	= 1e-1*dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;
	fcheckFrc 	= fcheck;


	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);

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
			cout << "		M E S O P H Y L L 					" << endl;
			cout << " 	F I R E 								" << endl;
			cout << "		E N T H A L P Y 					" << endl;
			cout << "	M I N I M I Z A T I O N 				" << endl;
			cout << "		W I T H   							" << endl;
			cout << "   L E E S  E D W A R D S      			" << endl;
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
			cout << "	** gamma 	= " << gamma << endl;
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
				x[i] -= L[i % NDIM]*floor(x[i]/L[i % NDIM]);
				if (i % NDIM == 0)
					x[i] -= floor(x[i+1]/L[1])*gamma*L[0];

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
		for (i=0; i<vertDOF; i++){
			// update position
			x[i] += dt*(v[i] + 0.5*x[i]*(Pi/V));

			// box replacement because of Lees Edwards
			x[i] -= L[i % NDIM]*floor(x[i]/L[i % NDIM]);

			// Lees-Edwards
			if (i % NDIM == 0)
				x[i] -= floor(x[i+1]/L[1])*gamma*L[0];
		}
		V += dt*Pi;
		L[0] = sqrt(V);
		L[1] = sqrt(V);
		for (d = 0; d < NDIM; d++)
			lb[d] = L[d] / sb[d];

		// update forces at fixed shear strain, FIXED CONTACT NETWORK
		mesoNetworkForceUpdate(gamma, gijtmp);

		// update instantaneous pressure
		P = 0.5*(stress[0] + stress[1]);

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
		cout << "		M E S O P H Y L L 					" << endl;
		cout << " 	F I R E 								" << endl;
		cout << "		E N T H A L P Y 					" << endl;
		cout << "	M I N I M I Z A T I O N 				" << endl;
		cout << "		W I T H   							" << endl;
		cout << "   L E E S  E D W A R D S      			" << endl;
		cout << "		C O N V E R G E D! 					" << endl;
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
		cout << "	** gamma 	= " << gamma << endl;
		cout << endl << endl;
	}
}


// use FIRE to minimize potential energy in mesophyll cells sinking toward center of box
void meso2D::mesoPinFIRE(vector<double> &xpin, double Ftol, double dt0, double kcspring){
	// local variables
	int i, ci, vi, gi, k;
	double cxtmp, cytmp;

	// check to see if cell linked-list has been initialized
	if (NBX == -1){
		cerr << "	** ERROR: In meso2D::mesoPinFIRE, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
		exit(1);
	}

	// also check that xpin has been created
	if (xpin.size() < 1){
		cerr << "	** ERROR: in meso2D::mesoPinFIRE, xpin size = " << xpin.size() << ", has not yet been created. Ending. " << endl;
		exit(1);
	}


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

	dtmax   	= 20.0*dt;
	dtmin   	= 1e-1*dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;

	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);

	// relax forces using FIRE
	while ((fcheck > Ftol || fireit < NDELAY) && fireit < itmax){
		// compute P
		P = 0.0;
		for (i=0; i<vertDOF; i++)
			P += v[i]*F[i];

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "		M E S O P H Y L L 			" << endl;
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
		for (i=0; i<vertDOF; i++)
			x[i] += dt*v[i];

		// update forces in mesophyll network 
		mesoPinForceUpdate(xpin, kcspring);

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
		cout << "		M E S O P H Y L L 			" << endl;
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


// mesophyll network in the mesophyll ensemble
void meso2D::mesoNetworkNVE(ofstream &enout, meso2DMemFn forceCall, double T, double dt0, int NT, int NPRINTSKIP){
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
		for (i=0; i<vertDOF; i++)
			x[i] += dt*v[i];

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
			cout << "===========================================" << endl;
			cout << "		M E S O P H Y L L 			" << endl;
			cout << " 			 						" << endl;
			cout << "		N V E 						" << endl;
			cout << "===========================================" << endl;
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


// NVE at fixed shear strain
void meso2D::mesoShearStrainNVE(ofstream &enout, double gamma, double T, double dt0, int NT, int NPRINTSKIP, vector<bool> &gijtmp){
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
		for (i=0; i<vertDOF; i++)
			x[i] += dt*v[i];

		// FORCE UPDATE
		mesoNetworkForceUpdate(gamma, gijtmp);

		// compute kinetic energy
		K = vertexKineticEnergy();

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*F[i]*dt;

		// update sim clock
		simclock += dt;

		// print to console and file
		if (t % NPRINTSKIP == 0){

			// print to console
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "		M E S O P H Y L L 			" << endl;
			cout << " 			 						" << endl;
			cout << "		N V E 						" << endl;
			cout << "===========================================" << endl;
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
				printMesoShearConfigCTCS2D(gamma);
		}
	}
}








/******************************

	M E S O P H Y L L

	P R O T O C O L S

*******************************/


// simulate network extension using decompression
void meso2D::mesoNetworkExtension(meso2DMemFn forceCall, double Ftol, double dt0, double delShrink, double dphiPrint, double phiMin){
	// local variables
	int ci, k = 0;
	double lastPrintPhi, scaleFactor = 1.0 - 2.0*delShrink;

	// current preferred packing fraction (neglect contribution from vertices)
	double phi0 = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		phi0 += a0.at(ci);
	phi0 /= L[0]*L[1];
	lastPrintPhi = phi0;

	// vector of vertices at edge of interfaces
	vector<bool> edge_verts;

	// loop until phi0 < phiMin
	while (phi0 > phiMin && k < itmax){
		// relax current configuration
		mesoFIRE(forceCall, Ftol, dt0);

		// find vertices at edge of interfaces
		edge_verts.resize(NVTOT);
		fill(edge_verts.begin(),edge_verts.end(),0);
		findInterfaceEdges(edge_verts);

		// break contact network
		updateMesophyllBondNetwork(edge_verts);

		// age particle shapes
		dt = delShrink;
		ageMesophyllShapeParameters();

		// add vertices
		if (NVTOT < NVMAX)
			addMesophyllCellMaterial(cL);


		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << "												" << endl;
		cout << " 	M E S O P H Y L L 							" << endl;
		cout << "	 											" << endl;
		cout << " 	N E T W O R K   							" << endl;
		cout << "												" << endl;
		cout << "	E X T E N S I O N 							" << endl;
		cout << " 												" << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi0 			= " << phi0 << endl;
		cout << "	* phi 			= " << vertexPackingFraction2D() << endl;
		cout << "	* phi0 w/ verts = " << vertexPreferredPackingFraction2D() << endl;
		cout << "	* lastPrintPhi 	= " << lastPrintPhi << endl;
		cout << "	* P 			= " << Pinst << endl;
		cout << "	* U 		 	= " << U << endl << endl;
		cout << "	* Contacts:" << endl;
		cout << "	* Nvv 			= " << vvContacts() << endl;
		cout << "	* Ncc 			= " << ccContacts() << endl << endl;
		cout << "	* Aging:" << endl;
		cout << "	* meanl0 		= " << meanl0() << endl;
		cout << "	* meancalA0 	= " << meancalA0() << endl;
		cout << "	* meant0 		= " << meant0() << endl;
		cout << "	* meankb 		= " << meankb() << endl;
		cout << endl << endl;

		// print positions if change in packing fraction is large enough
		if ((lastPrintPhi - phi0) > dphiPrint){
			// print positions
			printMesoNetwork2D();

			// update last phi when printed, for next time
			lastPrintPhi = phi0;
		}

		// scale particle sizes
		scaleParticleSizes2D(scaleFactor);

		// update new preferred fraction
		phi0 = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			phi0 += a0.at(ci);
		phi0 /= L[0]*L[1];

		// update iterate
		k++;
	}
}

// drag cell pins away from center (only cells with indices > cellskip)
void meso2D::mesoPinExtension(double Ftol, double dt0, double hmax, double dh, double dhprint, double kcspring, int cellskip){
	// local variables
	int k=0, ci, gi, vi;
	double phi0, cx, cy, dcx, dcy, h=0.0, lastPrinth=-10.0;
	vector<double> th(NCELLS,0.0);
	vector<double> xpin(NDIM*NCELLS,0.0);

	// determine extension direction angles + initial pin placement
	for (ci=0; ci<NCELLS; ci++){
		// get center of mass
		com2D(ci,cx,cy);

		// get distances to box center
		dcx = cx - 0.5*L[0];
		if (pbc[0])
			dcx -= L[0]*round(dcx/L[0]);

		dcy = cy - 0.5*L[1];
		if (pbc[1])
			dcy -= L[1]*round(dcy/L[1]);

		// get angle
		th[ci] = atan2(dcy,dcx);

		// save pin location
		xpin[NDIM*ci] = cx;
		xpin[NDIM*ci + 1] = cy;
	}

	// vector of vertices at edge of interfaces
	vector<bool> edge_verts;


	// loop over extension steps
	while (h < hmax && k < itmax){
		// printMesoPin2D(xpin, h);
		// printMesoBondNetwork();
		// if (k > 5)
		// 	exit(1);

		// relax current configuration
		setdt(dt0);
		mesoPinFIRE(xpin,Ftol,dt0,kcspring);

		// find vertices at edge of interfaces
		edge_verts.resize(NVTOT);
		fill(edge_verts.begin(),edge_verts.end(),0);
		findInterfaceEdges(edge_verts);

		// update contact network
		updateMesophyllBondNetwork(edge_verts);

		// age particle shapes
		dt = dh;
		ageMesophyllShapeParameters();

		// add vertices
		if (NVTOT < NVMAX)
			addMesophyllCellMaterial(cL);

		
		// relaxByAdding();

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << "												" << endl;
		cout << " 	M E S O P H Y L L 							" << endl;
		cout << "	 											" << endl;
		cout << " 	N E T W O R K . P I N   					" << endl;
		cout << "												" << endl;
		cout << "	E X T E N S I O N 							" << endl;
		cout << " 												" << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* NVTOT / NVMAX		= " << NVTOT << "/" << NVMAX << endl << endl;
		cout << "	* phi0 			= " << phi0 << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* h 			= " << h << endl;
		cout << "	* lastPrinth 	= " << lastPrinth << endl;
		cout << "	* P 			= " << 0.5*(stress[0] + stress[1]) << endl;
		cout << "	* U 		 	= " << U << endl << endl;
		cout << "	* Contacts:" << endl;
		cout << "	* Nvv 			= " << vvContacts() << endl;
		cout << "	* Ncc 			= " << ccContacts() << endl << endl;
		cout << "	* Aging:" << endl;
		cout << "	* meanl0 		= " << meanl0() << endl;
		cout << "	* meancalA0 	= " << meancalA0() << endl;
		cout << "	* meant0 		= " << meant0() << endl;
		cout << "	* meankb 		= " << meankb() << endl;
		cout << endl << endl;

		// print positions if change in packing fraction is large enough
		if (abs(lastPrinth - h) > dhprint){
			// print positions and bonds
			printMesoPin2D(xpin, h);
			printMesoBondNetwork();

			// update last print h
			lastPrinth = h;
		}

		// move pins + affine displace
		for (ci=0; ci<NCELLS; ci++){
			if (ci > cellskip){
				xpin[NDIM*ci] += dh*cos(th[ci]);
				xpin[NDIM*ci + 1] += dh*sin(th[ci]);

				gi = szList[ci];
				for (vi=0; vi<nv[ci]; vi++){
					x[NDIM*gi] += dh*cos(th[ci]);
					x[NDIM*gi + 1] += dh*sin(th[ci]);
					gi = gi + 1;
				}
			}
		}

		// scale particle sizes
		// scaleParticleSizes2D(1.0 - dh);

		// change particle sizes to relax pressure
		

		// update new packing fraction
		phi0 = vertexPreferredPackingFraction2D();

		// update new h
		h += dh;

		// update iterate
		k++;
	}
}

// grow network without boundary
void meso2D::mesoFreeGrowth(meso2DMemFn forceCall, double Ftol, double dt0, double dl0, double dphi0, double dphiPrint, double a0max){
	// local variables
	int ci, gi, k = 0;
	double lastPrintPhi, scaleFactor;

	// current preferred packing fraction (neglect contribution from vertices)
	double phi0 = 0.0;
	for (ci=0; ci<NCELLS; ci++)
		phi0 += a0.at(ci);
	phi0 /= L[0]*L[1];
	lastPrintPhi = phi0;
	printMesoNetwork2D();
	printMesoBondNetwork();

	// determine max areas
	double a0check = a0.at(0);
	double a0lim = a0max*a0check;

	// edge vertices
	vector<bool> edge_verts;

	// loop over growth
	while (a0check < a0lim && k < itmax){
		// relax current configuration
		mesoFIRE(forceCall, Ftol, dt0);

		// find vertices at edge of interfaces
		edge_verts.resize(NVTOT);
		fill(edge_verts.begin(),edge_verts.end(),0);
		findInterfaceEdges(edge_verts);

		// break contact network
		updateMesophyllBondNetwork(edge_verts);

		// age particle shapes
		dt = dphi0;
		ageMesophyllShapeParameters();

		// add vertices
		if (NVTOT < NVMAX)
			addMesophyllCellMaterial(dl0);

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << "												" << endl;
		cout << " 	M E S O P H Y L L 							" << endl;
		cout << "	 											" << endl;
		cout << " 	F R E E   									" << endl;
		cout << "												" << endl;
		cout << "	G R O W T H 								" << endl;
		cout << " 												" << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi0 			= " << phi0 << endl;
		cout << " 	* a0check 		= " << a0check << endl;
		cout << " 	* lastPrintPhi 	= " << lastPrintPhi << endl;
		cout << "	* P 			= " << 0.5*(stress[0] + stress[1]) << endl;
		cout << "	* U 		 	= " << U << endl << endl;
		cout << "	* Contacts:" << endl;
		cout << "	* Nvv 			= " << vvContacts() << endl;
		cout << "	* Ncc 			= " << ccContacts() << endl << endl;
		cout << "	* Aging:" << endl;
		cout << "	* meanl0 		= " << meanl0() << endl;
		cout << "	* meancalA0 	= " << meancalA0() << endl;
		cout << "	* meant0 		= " << meant0() << endl;
		cout << "	* meankb 		= " << meankb() << endl;
		cout << endl << endl;

		// print positions if change in packing fraction is large enough
		if (abs(lastPrintPhi - phi0) > dphiPrint){
			// print positions
			printMesoNetwork2D();
			printMesoBondNetwork();

			// update last phi when printed, for next time
			lastPrintPhi = phi0;
		}

		// increase perimeters near void
		for (gi=0; gi<NVTOT; gi++){
			if (zv[gi] == 0 || zv[ip1[gi]] == 0)
				l0[gi] += (dl0*dphi0/NCELLS)*l0[gi];
		}


		// increase areas, update check
		phi0 = 0.0;
		for (ci=0; ci<NCELLS; ci++)
			phi0 += a0.at(ci);
		phi0 /= L[0]*L[1];

		// determine scale factor
		scaleFactor = sqrt((phi0 + dphi0)/phi0);

		// scale particle sizes
		scaleParticleSizes2D(scaleFactor);
		a0check = a0[0];


		// update iterate
		k++;
	}
}

// simulate network *formation* using *enthalpy minimization*
// Updates (02/14/22)
// 	* Changed perimeter growth to be on l0_im1, l0i only if vertex zv=0
// 	* Changed angle growth to be driven toward negative minimum value
// 
// 	* ISSUE: how to get homogenous growth of perimeter? Seems like 
// 	some regions "run off" while others stay small, causes heterogenous growth dynamics
// 	-- Could be fixed by growth of void bubbles
// 	-- SEEMS TO BE DRIVEN by non-homogenous contact breaking ... make it more so by labelling edge of voids differently? MAKE contacts break between different interfaces
// 	-- 02/15: Getting closer, but will get to circle stage and then poro starts decreasing ... need to make sure that in circle stage, perimeter growth drives boundary growth?
// 		** Try plotting which vertices are labelled as "edge", and 
void meso2D::mesoNetworkEnthalpyMin(meso2DMemFn forceCall, double Ftol, double dPtol, double dt0, double da0, double dl0, double t0_min, double P0, double phiMin, int NMINSKIP){
	// local variables
	int i, ci, gi, vi, d, k = 0;
	double G, a0old, phi, phi0, calAtmp;

	// initialize packing fraction definitions
	phi0 = vertexPreferredPackingFraction2D();
	phi = vertexPackingFraction2D();

	// list of edge contacts
	vector<bool> edge_verts;

	// find vertices at edge of contacts
	edge_verts.resize(NVTOT);
	fill(edge_verts.begin(),edge_verts.end(),0);
	findInterfaceEdges(edge_verts);

	// loop until phi0 < phiMin
	while (phi > phiMin && k < itmax){

		// // growth driven by ALL VOID VERTICES
		// for (gi=0; gi<NVTOT; gi++){
		// 	// global radii growth
		// 	r[gi] *= 1 + dl0*da0;

		// 	// drive void-facing perimeter growth
		// 	if (zv[gi] <= 0){
		// 		// grow segment im1 and i if vertex i not in contact
		// 		l0[im1[gi]] *= (1 + dl0);
		// 		l0[gi] *= (1.0 + dl0);

		// 		// drive curvature towards t0_min
		// 		if (t0[gi] > t0_min + dl0*cB)
		// 			t0[gi] -= dl0*cB;
		// 	}
		// 	// otherwise, relax angle toward 0 (flat interfaces)
		// 	else
		// 		t0[gi] *= 1 - dl0*cB;
		// }

		// // growth driven by EDGE VERTICES
		// for (gi=0; gi<NVTOT; gi++){
		// 	// global radii growth
		// 	r[gi] *= 1 + dl0*da0;

		// 	// growth between two edge vertices
		// 	if ( (edge_verts[gi] && (edge_verts[ip1[gi]] || zv[ip1[gi]] <= 0)) || (zv[gi] <= 0 && edge_verts[ip1[gi]]) ) {
		// 		// grow segment between edge vertices
		// 		l0[gi] *= (1.0 + dl0);

		// 		// drive curvature towards t0_min
		// 		if (t0[gi] > t0_min + dl0*cB)
		// 		// if (t0[gi] > -0.5*PI)
		// 			t0[gi] -= dl0*cB;

		// 		if (t0[ip1[gi]] > t0_min + dl0*cB)
		// 		// if (t0[gi] > -0.5*PI)
		// 			t0[ip1[gi]] -= dl0*cB;
		// 	}
		// 	else if (zv[gi] <= 0){
		// 		if (t0[gi] > t0_min + dl0*cB)
		// 			t0[gi] -= dl0*cB;
		// 	}
		// 	else
		// 		t0[gi] *= 1 - dl0*cB;
		// }

		// // area-level growth
		// for (ci=0; ci<NCELLS; ci++)
		// 	a0[ci] *= pow(1 + dl0*da0,2.0);


		// age and grow
		ageMesophyllShapeParameters(edge_verts, dl0, da0, t0_min);

		// relax current configuration
		mesoEnthalpyFIRE(forceCall, Ftol, dPtol, P0, dt0);

		// for (ci=0; ci<NCELLS; ci++)
		// 	a0[ci] = area(ci);
		// t0ToCurrent();

		// print MINIMIZED positions if change in packing fraction is large enough
		if (k % NMINSKIP == 0){
			printMesoNetworkCTCS2D();
			if (hessout.is_open())
				mesoPrintLinearResponse(forceCall, Ftol, dt0, P0);
		}

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << "												" << endl;
		cout << " 	M E S O P H Y L L 							" << endl;
		cout << "	 											" << endl;
		cout << " 	N E T W O R K   							" << endl;
		cout << "												" << endl;
		cout << "	E X T E N S I O N 							" << endl;
		cout << " 												" << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi0 			= " << phi0 << endl;
		cout << "	* phi 			= " << vertexPackingFraction2D() << endl;
		cout << "	* phi0 w/ verts = " << vertexPreferredPackingFraction2D() << endl;
		cout << "	* P 			= " << Pinst << endl;
		cout << "	* U 		 	= " << U << endl << endl;
		cout << "	* Contacts:" << endl;
		cout << "	* Nvv 			= " << vvContacts() << endl;
		cout << "	* Ncc 			= " << ccContacts() << endl << endl;
		cout << "	* Aging:" << endl;
		cout << "	* meanl0 		= " << meanl0() << endl;
		cout << "	* meancalA0 	= " << meancalA0() << endl;
		cout << "	* meant0 		= " << meant0() << endl;
		cout << "	* meankb 		= " << meankb() << endl;
		cout << endl << endl;

		// break contact network
		updateMesophyllBondNetwork(edge_verts);

		// recompute bond information
		// computeZ();

		// find vertices at edge of contacts
		edge_verts.resize(NVTOT);
		fill(edge_verts.begin(),edge_verts.end(),0);
		findInterfaceEdges(edge_verts);

		// // age particle shapes
		// dt = dl0;
		// ageMesophyllShapeParameters();

		// add vertices
		if (NVTOT < NVMAX)
			addMesophyllCellMaterial(edge_verts);

		// update packing fraction
		phi0 = vertexPreferredPackingFraction2D();
		phi = vertexPackingFraction2D();

		// update iterate
		k++;
	}
}




// drive voids by mesophyll void bubble growth
void meso2D::mesoBubbleEnthalpyMin(meso2DMemFn forceCall, double Ftol, double dPtol, double dt0, double da0, double P0, double phiMin, int NMINSKIP){
	// local variables
	bool noBubbles = 1;
	int NTOT, NTCELLS, i, ci, gi, vi, d, k = 0;
	double G, a0old, phi, phi0, calAtmp;
	double bubblePhi, growPhi;

	// list of edge contacts
	vector<bool> edge_verts;

	// find vertices at edge of contacts
	edge_verts.resize(NVTOT);
	fill(edge_verts.begin(),edge_verts.end(),0);
	findInterfaceEdges(edge_verts);

	// change number of cells, to be switched when voids are large enough
	NTOT = NCELLS;
	NTCELLS = NTOT - NBUBBLES;
	bubblePhi = 0.0;
	for (ci=NTCELLS; ci<NCELLS; ci++)
		bubblePhi += area(ci)/(L[0]*L[1]);
	NCELLS = NTCELLS;
	phi = vertexPackingFraction2D();
	// growPhi = phi - 2.0*bubblePhi;
	growPhi = 0.0;

	// initialize packing fraction definitions
	phi0 = vertexPreferredPackingFraction2D();
	phi = vertexPackingFraction2D();

	// loop until phi0 < phiMin
	while (phi > phiMin && k < 5){

		// relax current configuration
		mesoEnthalpyFIRE(forceCall, Ftol, dPtol, P0, dt0);

		// void bubbles
		if (noBubbles){
			// relax positions in fixed geometry
			relaxSmallBubbles();

			// change if phi is right
			if (phi < growPhi){
				noBubbles = 0;
				NCELLS = NTOT;
			}
		}
		else{
			// homogenously grow cells
			for (ci=NTCELLS; ci<NCELLS; ci++){
				a0[ci] *= pow(1 + da0,2.0);
				for (vi=0; vi<nv[ci]; vi++){
					gi = szList[ci] + vi;
					l0[gi] *= 1 + da0;
					r[gi] *= 1 + da0;
				}
			}
		}

		// change to cells
		for (ci=0; ci<NTCELLS; ci++)
			a0[ci] *= pow(1 + da0,2.0);

		// print MINIMIZED positions if change in packing fraction is large enough
		if (k % NMINSKIP == 0){
			printMesoNetworkCTCS2D();
			if (hessout.is_open())
				mesoPrintLinearResponse(forceCall, Ftol, dt0, P0);
		}

		// output to console
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << "===============================================" << endl << endl;
		cout << "												" << endl;
		cout << " 	M E S O P H Y L L 							" << endl;
		cout << "	 											" << endl;
		cout << " 	N E T W O R K   							" << endl;
		cout << "												" << endl;
		cout << "	E X T E N S I O N 							" << endl;
		cout << " 												" << endl;
		cout << "===============================================" << endl;
		cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
		cout << endl;
		cout << "	* k 			= " << k << endl;
		cout << "	* phi0 			= " << phi0 << endl;
		cout << "	* phi 			= " << vertexPackingFraction2D() << endl;
		cout << "	* phi0 w/ verts = " << vertexPreferredPackingFraction2D() << endl;
		cout << "	* P 			= " << Pinst << endl;
		cout << "	* U 		 	= " << U << endl << endl;
		cout << "	* Contacts:" << endl;
		cout << "	* Nvv 			= " << vvContacts() << endl;
		cout << "	* Ncc 			= " << ccContacts() << endl << endl;
		cout << "	* Aging:" << endl;
		cout << "	* meanl0 		= " << meanl0() << endl;
		cout << "	* meancalA0 	= " << meancalA0() << endl;
		cout << "	* meant0 		= " << meant0() << endl;
		cout << "	* meankb 		= " << meankb() << endl;
		cout << endl << endl;

		// relax shape parameters in presence of void bubbles
		for (ci=0; ci<NCELLS; ci++){
			for (vi=0; vi<nv[ci]; vi++){
				gi = szList[ci] + vi;
				l0[gi] = perimeter(ci)/nv[ci];
			}
		}

		// break contact network
		updateMesophyllBondNetwork(edge_verts);

		// find vertices at edge of contacts
		edge_verts.resize(NVTOT);
		fill(edge_verts.begin(),edge_verts.end(),0);
		findInterfaceEdges(edge_verts);

		// add vertices
		if (NVTOT < NVMAX)
			addMesophyllCellMaterial(edge_verts);

		// update packing fraction
		phi0 = vertexPreferredPackingFraction2D();
		phi = vertexPackingFraction2D();

		// update iterate
		k++;
	}
}




/******************************

	M E S O P H Y L L

	P R O T O C O L

	H E L P E R S

*******************************/


// update mesophyll network using MC
void meso2D::updateMesophyllBondNetwork(vector<bool> &edge_verts){
	// local variables
	bool isConnected, canBreak;
	int imove, jmove, gitmp;
	int cijctc, zitmp, zjtmp, ci, cj, ck, vi, vj, gi, gj, hi, hj;
	double dx, dy, sij, rij, zij, dU, poff, h=ctch, h2=h*h, rdraw;
	double lix, liy, lim1x, lim1y, ljx, ljy, ljm1x, ljm1y, btmp;

	// loop over pairs of vertices, check whether to connect or detach
	for (ci=0; ci<NCELLS; ci++){
		// get total number of cells in bonded contact with ci
		zitmp = mesoBondedCTCS(ci);

		// loop over other cells to determine which bonds to break
		for (cj=(ci+1); cj<NCELLS; cj++){
			// get total number of cells in bonded contact with cj, also count bonded contacts between ci and cj
			zjtmp = mesoBondedCTCS(cj);
			cijctc = mesoBondedPAIRS(ci,cj);

			// zij: determines strength of bond attraction
			zij = 0.5*(zc[ci] + zc[cj])*ctcdel + 1.0;

			// use global label of vi
			gi = szList[ci];
			for (vi=0; vi<nv[ci]; vi++){
				// use global label of vj
				gj = szList[cj];
				for (vj=0; vj<nv[cj]; vj++){
					// check if connected, otherwise go to next vertex
					isConnected = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];

					if (isConnected){
						// contact distance
						sij = r[gi] + r[gj];

						// get vertex-vertex distances
						dx = x[NDIM*gj] - x[NDIM*gi];
						if (pbc[0])
							dx -= L[0]*round(dx/L[0]);

						dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
						if (pbc[1])
							dy -= L[1]*round(dy/L[1]);

						// true distance
						rij = sqrt(dx*dx + dy*dy);

						// check if bond should move, if so don't break
						imove = 0;
						jmove = 0;



						// bond vectors in x-direction
						lim1x = x[NDIM*gi] - x[NDIM*im1[gi]];
						lix = x[NDIM*ip1[gi]] - x[NDIM*gi];
						ljm1x = x[NDIM*gj] - x[NDIM*im1[gj]];
						ljx = x[NDIM*ip1[gj]] - x[NDIM*gj];
						if (pbc[0]){
							lim1x -= L[0]*round(lim1x/L[0]);
							lix -= L[0]*round(lix/L[0]);
							ljm1x -= L[0]*round(ljm1x/L[0]);
							ljx -= L[0]*round(ljx/L[0]);
						}

						// bond vectors in y-direction
						lim1y = x[NDIM*gi + 1] - x[NDIM*im1[gi] + 1];
						liy = x[NDIM*ip1[gi] + 1] - x[NDIM*gi + 1];
						ljm1y = x[NDIM*gj + 1] - x[NDIM*im1[gj] + 1];
						ljy = x[NDIM*ip1[gj] + 1] - x[NDIM*gj + 1];
						if (pbc[0]){
							lim1y -= L[0]*round(lim1y/L[0]);
							liy -= L[0]*round(liy/L[0]);
							ljm1y -= L[0]*round(ljm1y/L[0]);
							ljy -= L[0]*round(ljy/L[0]);
						}



						// // -- check whether or not to move

						// // new bond length moving from i to ip1
						// btmp = pow(dx - lix,2.0) + pow(dy - liy,2.0);
						// if (btmp < rij*rij && !gij[NVTOT*ip1[gi] + gj - (ip1[gi]+1)*(ip1[gi]+2)/2]){
						// 	// set has moved to +1
						// 	imove = 1;

						// 	// switch v-v bond to gi+1
						// 	gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 0;
						// 	gij[NVTOT*ip1[gi] + gj - (ip1[gi]+1)*(ip1[gi]+2)/2] = 1;

						// 	// update vertex-vertex distances
						// 	dx = x[NDIM*gj] - x[NDIM*ip1[gi]];
						// 	if (pbc[0])
						// 		dx -= L[0]*round(dx/L[0]);

						// 	dy = x[NDIM*gj + 1] - x[NDIM*ip1[gi] + 1];
						// 	if (pbc[1])
						// 		dy -= L[1]*round(dy/L[1]);

						// 	// true distance
						// 	rij = sqrt(dx*dx + dy*dy);
						// }

						// // moving from i to im1
						// btmp = pow(dx + lim1x,2.0) + pow(dy + lim1y,2.0);
						// if (btmp < rij*rij && imove == 0 && !gij[NVTOT*im1[gi] + gj - (im1[gi]+1)*(im1[gi]+2)/2]){
						// 	// set has moved to -1
						// 	imove = -1;

						// 	// switch v-v bond to gi-1
						// 	gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 0;
						// 	gij[NVTOT*im1[gi] + gj - (im1[gi]+1)*(im1[gi]+2)/2] = 1;

						// 	// update vertex-vertex distances
						// 	dx = x[NDIM*gj] - x[NDIM*im1[gi]];
						// 	if (pbc[0])
						// 		dx -= L[0]*round(dx/L[0]);

						// 	dy = x[NDIM*gj + 1] - x[NDIM*im1[gi] + 1];
						// 	if (pbc[1])
						// 		dy -= L[1]*round(dy/L[1]);

						// 	// true distance
						// 	rij = sqrt(dx*dx + dy*dy);
						// }


						// ALSO CHECK FOR MOVING J TO NEIGHBORING VERTICES!!

						// // moving from j to jp1
						// btmp = pow(dx + ljx,2.0) + pow(dy + ljy,2.0);
						// if (btmp < rij*rij){
						// 	// set jmove = 1
						// 	jmove = 1;

						// 	// determine if need to use new gi
						// 	if (imove == 1)
						// 		gitmp = ip1[gi];
						// 	else if (imove == -1)
						// 		gitmp = im1[gi];
						// 	else
						// 		gitmp = gi;

						// 	// switch v-v bond to gitmp & ip1[gj]
						// 	gij[NVTOT*gitmp + gj - (gitmp+1)*(gitmp+2)/2] = 0;
						// 	gij[NVTOT*gitmp + ip1[gj] - (gitmp+1)*(gitmp+2)/2] = 1;
						// }

						// // moving from j to jm1
						// btmp = pow(dx - ljm1x,2.0) + pow(dy - ljm1y,2.0);
						// if (btmp < rij*rij && jmove == 0){
						// 	// set jmove to -1
						// 	jmove = -1;

						// 	// determine if need to use new gi
						// 	if (imove == 1)
						// 		gitmp = ip1[gi];
						// 	else if (imove == -1)
						// 		gitmp = im1[gi];
						// 	else
						// 		gitmp = gi;

						// 	// switch v-v bond to gitmp & ip1[gj]
						// 	gij[NVTOT*gitmp + gj - (gitmp+1)*(gitmp+2)/2] = 0;
						// 	gij[NVTOT*gitmp + im1[gj] - (gitmp+1)*(gitmp+2)/2] = 1;
						// }

						// only check bond if extended & has not moved
						if (rij > sij){
							// change in energy from bond breaking
							dU = 1.0 - 0.5*(pow(1 - (rij/sij),2.0)/h2);

							// only break stochastically at edges
							canBreak = (edge_verts[gi] && edge_verts[gj]);
							// canBreak = 1;

							// remove if bond detaching decreases energy
							if (dU < 0){
								// remove vertex-vertex bond
								gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 0;

								// remove single cell-cell contact
								cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]--;

								// update contact lists 
								zc[ci]--;
								zc[cj]--;

								zv[gi]--;
								zv[gj]--;
							}
							else{
								// else, remove conditionally
								poff = exp(-(betaEff*dU)/zij);
								rdraw = drand48();

								// detach
								if (poff > rdraw && canBreak){
									// detach vv contact
									gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 0;

									// remove single cell-cell contact
									cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]--;

									// update contact lists 
									zc[ci]--;
									zc[cj]--;

									zv[gi]--;
									zv[gj]--;
								}
							}
						}
					}
					// increment global vertex index
					gj++;
				}
				

				// increment global vertex index
				gi++;
			}
		}
	}
}

// age mesophyll shape parameters
void meso2D::ageMesophyllShapeParameters(vector<bool> &edge_verts, double dl0, double da0, double t0_min){
	// local variables
	int gi, ci, vi;
	double lix, liy, lim1x, lim1y, li, ti, sini, cosi;
	double dl0_tmp;

	// count number of void vertices per cell
	vector<int> nvoid(NCELLS,0);

	// grow areas, radii
	gi = 0; 
	for (ci=0; ci<NCELLS; ci++){
		// grow area
		a0[ci] *= pow(1 + dl0*da0,2.0);
		for (vi=0; vi<nv[ci]; vi++){
			// grow vertex radius
			r[gi] *= (1 + da0*dl0);

			// count number of void segments
			if (zv[gi] == 0)
				nvoid[ci]++;

			// increment global vertex counter
			gi++;
		}
	}

	// grow and age perimeter
	for (gi=0; gi<NVTOT; gi++){
		// cell index
		cindices(ci,vi,gi);
		dl0_tmp = dl0*(1.0/(1.0 + (nvoid[ci]/nv[ci])));

		// segment from i to ip1
		lix = x[NDIM*ip1[gi]] - x[NDIM*gi];
		if (pbc[0])
			lix -= L[0]*round(lix/L[0]);

		liy = x[NDIM*ip1[gi] + 1] - x[NDIM*gi + 1];
		if (pbc[1])
			liy -= L[1]*round(liy/L[1]);

		// segment from im1 to i
		lim1x = x[NDIM*gi] - x[NDIM*im1[gi]];
		if (pbc[0])
			lim1x -= L[0]*round(lim1x/L[0]);

		lim1y = x[NDIM*gi + 1] - x[NDIM*im1[gi] + 1];
		if (pbc[1])
			lim1y -= L[1]*round(lim1y/L[1]);

		// segment length
		li = sqrt(lix*lix + liy*liy);

		// local bending
		sini = lim1x*liy - lim1y*lix;
		cosi = lim1x*lix + lim1y*liy;
		ti = atan2(sini,cosi);

		// // age near edges
		// if ( (edge_verts[gi] && (edge_verts[ip1[gi]] || zv[ip1[gi]] == 0)) || (zv[gi] == 0 && edge_verts[ip1[gi]]) ){
		// 	l0[gi] *= (1.0 + dl0);
		// 	if (t0[gi] > t0_min + dl0*cB)
		// 		t0[gi] -= dl0*cB;
		// }
		// // grow on void
		// else if (zv[gi] == 0){
		// 	l0[gi] += dl0*cL*(li - l0[gi]);
		// 	t0[gi] += dl0*cB*(ti - t0[gi]);
		// }

		// grow along void areas
		if (zv[gi] == 0){
			l0[gi] *= (1.0 + dl0_tmp);
			if (t0[gi] > t0_min + dl0_tmp*cB)
				t0[gi] -= dl0*cB;
		}
		// if ctc, age toward 0
		else{
			t0[gi] *= 1 - dl0*cB;
			l0[gi] += dl0*(li - l0[gi]);
		}
	}
}


// add vertices between edge contacts or when distance between two void vertices is too large
void meso2D::addMesophyllCellMaterial(vector<bool> &edge_verts){
	// check that edge verts is correct size
	if (edge_verts.size() != NVTOT){
		cout << "ERROR: in addMesophyllCellMaterial, edge_verts.size != NVTOT. size=" << edge_verts.size() << ", NVTOT=" << NVTOT << "...ending. " << endl;
		exit(1);
	}

	// local variables
	bool gtmp = 0;
	int d, gi, gj, gim1, gip1, ctc_im_jn, ctc_ip1m_kp, ci, vi, cin, vin, cip1p, vip1p;
	double dx, dy, dli, dlim1, di, dmin, lix, liy, meanl, testl;

	// vector to pushback vertices where extra material is added
	vector<int> growthVerts;

	// loop over vertices, decide which to grow
	for (gi=0; gi<NVTOT; gi++){
		// info for vertex gi
		cindices(ci,vi,gi);
		gip1 = ip1[gi];
		gim1 = im1[gi];

		// check distance to backward neighbor
		dx = x[NDIM*gi] - x[NDIM*gim1];
		dy = x[NDIM*gi + 1] - x[NDIM*gim1 + 1];
		if (pbc[0])
			dx -= L[0]*round(dx/L[0]);
		if (pbc[1])
			dy -= L[1]*round(dy/L[1]);
		dlim1 = sqrt(dx*dx + dy*dy);

		// check distance to forward neighbor
		dx = x[NDIM*gip1] - x[NDIM*gi];
		dy = x[NDIM*gip1 + 1] - x[NDIM*gi + 1];
		if (pbc[0])
			dx -= L[0]*round(dx/L[0]);
		if (pbc[1])
			dy -= L[1]*round(dy/L[1]);
		dli = sqrt(dx*dx + dy*dy);

		// add if void vertex and over extended
		if (zv[gi] == 0){
			// if distance between two vertices is more than twice the vertex radius, add
			if (dlim1 > 4.0*r[gi])
				growthVerts.push_back(gim1);
			if (dli > 4.0*r[gi])
				growthVerts.push_back(gi);

		}
		// else if (edge_verts[gi] && edge_verts[gip1]){
		// 	growthVerts.push_back(gi);
		// }
		// add if neighbors connected to different cells, or if center not connected,
		// but ip1 and im1 are connected
		else if (zv[gi] > 0 && zv[gip1] > 0){
			// check min contact to gi
			dmin = 1e6;
			for (gj=0; gj<NVTOT; gj++){
				if (gj < gi)
					gtmp = gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2];
				else if (gj > gi)
					gtmp = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];

				if (gtmp && gi != gj){
					// distance between j and i in x
					dx = x[NDIM*gj] - x[NDIM*gi];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);

					// distance between j and i in y
					dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// full distance
					di = dx*dx + dy*dy;
					if (di < dmin){
						dmin = di;
						ctc_im_jn = gj;
					}
				}
			}

			// check min contact to gip1
			dmin = 1e6;
			for (gj=0; gj<NVTOT; gj++){
				if (gj < gip1)
					gtmp = gij[NVTOT*gj + gip1 - (gj+1)*(gj+2)/2];
				else if (gj > gip1)
					gtmp = gij[NVTOT*gip1 + gj - (gip1+1)*(gip1+2)/2]; 

				if (gtmp && gip1 != gj){
					// distance between j and i in x
					dx = x[NDIM*gj] - x[NDIM*gip1];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);

					// distance between j and i in y
					dy = x[NDIM*gj + 1] - x[NDIM*gip1 + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// full distance
					di = dx*dx + dy*dy;
					if (di < dmin){
						dmin = di;
						ctc_ip1m_kp = gj;
					}
				}
			}

			// if contacts are on different cells, add vertex between contacts if space available
			cindices(cin,vin,ctc_im_jn);
			cindices(cip1p,vip1p,ctc_ip1m_kp);
			if (cin != cip1p){
				cout << "gi=" << gi << ", cin=" << cin << ", gip1=" << gip1 << ", cip1p=" << cip1p << endl;
				growthVerts.push_back(gi);
			}
		}

		// add if lone vertex between two connected vertices, connected to different cells
		else if (zv[gim1] > 0 && zv[gip1] > 0 && zv[gi] == 0){
			// check min contact to gim1
			dmin = 1e6;
			for (gj=0; gj<NVTOT; gj++){
				if (gj < gim1)
					gtmp = gij[NVTOT*gj + gim1 - (gj+1)*(gj+2)/2];
				else if (gj > gim1)
					gtmp = gij[NVTOT*gim1 + gj - (gim1+1)*(gim1+2)/2];

				if (gtmp && gim1 != gj){
					// distance between j and i in x
					dx = x[NDIM*gj] - x[NDIM*gim1];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);

					// distance between j and i in y
					dy = x[NDIM*gj + 1] - x[NDIM*gim1 + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// full distance
					di = dx*dx + dy*dy;
					if (di < dmin){
						dmin = di;
						ctc_im_jn = gj;
					}
				}
			}

			// check min contact to gip1
			dmin = 1e6;
			for (gj=0; gj<NVTOT; gj++){
				if (gj < gip1)
					gtmp = gij[NVTOT*gj + gip1 - (gj+1)*(gj+2)/2];
				else if (gj > gip1)
					gtmp = gij[NVTOT*gip1 + gj - (gip1+1)*(gip1+2)/2]; 

				if (gtmp && gip1 != gj){
					// distance between j and i in x
					dx = x[NDIM*gj] - x[NDIM*gip1];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);

					// distance between j and i in y
					dy = x[NDIM*gj + 1] - x[NDIM*gip1 + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// full distance
					di = dx*dx + dy*dy;
					if (di < dmin){
						dmin = di;
						ctc_ip1m_kp = gj;
					}
				}
			}

			// if contacts are on different cells, add vertex between contacts
			cindices(cin,vin,ctc_im_jn);
			cindices(cip1p,vip1p,ctc_ip1m_kp);
			if (cin != cip1p){
				growthVerts.push_back(gim1);
				growthVerts.push_back(gi);
			}
		}
	}
	// Loop over growth locations, run growth protocol
	if (growthVerts.size() > 0){
		for (gi=0; gi<growthVerts.size(); gi++){
			cout << "adding vertex between gi=" << growthVerts.at(gi)+gi << " and gip1=" << ip1[growthVerts.at(gi)+gi] << endl;
			addVertex(growthVerts.at(gi)+gi,l0.at(growthVerts.at(gi)+gi));
		}
	}
	growthVerts.clear();
}


// find vertices at the edges of cell-cell interfaces
void meso2D::findInterfaceEdges(vector<bool> &edge_verts){
	// check that edge verts is correct size
	if (edge_verts.size() != NVTOT){
		cout << "ERROR: in findInterfaceEdges, edge_verts.size != NVTOT. size=" << edge_verts.size() << ", NVTOT=" << NVTOT << "...ending. " << endl;
		exit(1);
	}

	// local variables
	bool diff, gitmp, gjtmp;
	int gi, gim1, gip1, gj, gjp1, cthis, vthis, ci, cj, vi, vj;


	// loop over vertices, check for edges of cell-cell interfaces
	for (gi=0; gi<NVTOT; gi++){
		// only proceed if gi is in contact
		if (zv[gi] > 0){
			// neighboring vertices to check if on edge
			gip1 = ip1[gi];
			gim1 = im1[gi];

			// get cell info of gi
			cindices(cthis, vthis, gi);

			// check forward neighbor
			if (zv[gip1] == 0){
				edge_verts[gi] = 1;
				continue;
			}
			else{
				// find cells that gi and gip1 are in contact with
				diff = 0;
				for (ci=0; ci<NCELLS; ci++){
					if (ci != cthis){
						for (vi=0; vi<nv[ci]; vi++){
							// get possible contact with gi
							gj = szList[ci] + vi;

							// see if gi is connected to gj
							if (gj < gi)
								gitmp = gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2];
							else if (gj > gi)
								gitmp = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];

							// if gj a match, get possible contact with gip1
							if (gitmp){
								for (cj=0; cj<NCELLS; cj++){
									if (cj != cthis && cj != ci){
										for (vj=0; vj<nv[cj]; vj++){
											// get possible contact with gip1
											gjp1 = szList[cj] + vj;

											// see if gip1 is connected to gjp1
											if (gjp1 < gip1)
												gjtmp = gij[NVTOT*gjp1 + gip1 - (gjp1+1)*(gjp1+2)/2];
											else if (gjp1 > gip1)
												gjtmp = gij[NVTOT*gip1 + gjp1 - (gip1+1)*(gip1+2)/2];

											// if gip1 is in contact with gip1, on different cell than gj, break!
											if (gjtmp){
												diff = 1;
												// cout << "Found edge vertex! On (" << cthis << ", " << vthis << "): vi ctc to (" << ci << ", " << vi << ");   vip1 ctc to (" << cj << ", " << vj << ")..." << endl;
												break;
											}
										}
									}
									if (diff)
										break;
								}
							}
							if (diff)
								break;
						}
					}
					if (diff)
						break;
				}
				if (diff)
					edge_verts[gi] = 1;
			}

			// check backward neighbor (only need to check if not ctc...other case handled above)
			if (zv[gim1] == 0){
				edge_verts[gi] = 1;
				continue;
			}
		}
	}
}


// recompute number of vertices on each bond
void meso2D::computeZ(){
	// local variables
	int gi, ci, vi, gj, cj, vj;

	// reset z
	fill(zv.begin(),zv.end(),0);
	fill(zc.begin(),zc.end(),0);

	// loop over each bond, recompute zc and zv
	for (gi=0; gi<NVTOT; gi++){
		cindices(ci,vi,gi);
		for (gj=(gi+1); gj<NVTOT; gj++){
			cindices(cj,vj,gj);
			if (cj == ci)
				continue;

			// check if bonded
			if (gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2]){
				zv[gi]++;
				zv[gj]++;

				zc[ci]++;
				zc[cj]++;
			}
		}
	}
}


// age mesophyll shape parameters
void meso2D::ageMesophyllShapeParameters(){
	// local variables
	int gi;
	double lix, liy, lim1x, lim1y, li, ti, sini, cosi;

	for (gi=0; gi<NVTOT; gi++){	
		// segment from i to ip1
		lix = x[NDIM*ip1[gi]] - x[NDIM*gi];
		if (pbc[0])
			lix -= L[0]*round(lix/L[0]);

		liy = x[NDIM*ip1[gi] + 1] - x[NDIM*gi + 1];
		if (pbc[1])
			liy -= L[1]*round(liy/L[1]);

		// segment from im1 to i
		lim1x = x[NDIM*gi] - x[NDIM*im1[gi]];
		if (pbc[0])
			lim1x -= L[0]*round(lim1x/L[0]);

		lim1y = x[NDIM*gi + 1] - x[NDIM*im1[gi] + 1];
		if (pbc[1])
			lim1y -= L[1]*round(lim1y/L[1]);

		// -- age perimeter segments

		// segment length
		li = sqrt(lix*lix + liy*liy);

		// update l0 (either void or contact)
		if (l0[gi] < li){
			if (zv[gi] > 0 && zv[ip1[gi]] > 0)
				l0[gi] += (1.0-aL)*dt*cL*(li - l0[gi]);
			else
				l0[gi] += aL*dt*cL*(li - l0[gi]);
		}
		

		// -- age angles

		// angle trig functions
		sini = lim1x*liy - lim1y*lix;
		cosi = lim1x*lix + lim1y*liy;

		// angle
		ti = atan2(sini,cosi);

		// update t0
		if (zv[gi] > 0 && zv[ip1[gi]] > 0)
			t0[gi] -= dt*cB*t0[gi];
		else
			t0[gi] += dt*cB*(ti - t0[gi]);


		// age bending mechanical constant
		if (kbi[gi] < kbmax - cKb)
			kbi[gi] += dt*cKb;
	}
}


// relax via addition of material
void meso2D::relaxByAdding(){
	// local variables
	int gi, ngrow;
	double lix, liy, lim1x, lim1y, li, ti, sini, cosi;

	// vector to pushback vertices where extra material is added
	vector<int> growthVerts;

	ngrow = 0;
	for (gi=0; gi<NVTOT; gi++){	
		// segment from i to ip1
		lix = x[NDIM*ip1[gi]] - x[NDIM*gi];
		if (pbc[0])
			lix -= L[0]*round(lix/L[0]);

		liy = x[NDIM*ip1[gi] + 1] - x[NDIM*gi + 1];
		if (pbc[1])
			liy -= L[1]*round(liy/L[1]);

		// segment from im1 to i
		lim1x = x[NDIM*gi] - x[NDIM*im1[gi]];
		if (pbc[0])
			lim1x -= L[0]*round(lim1x/L[0]);

		lim1y = x[NDIM*gi + 1] - x[NDIM*im1[gi] + 1];
		if (pbc[1])
			lim1y -= L[1]*round(lim1y/L[1]);

		// -- age perimeter segments

		// segment length
		li = sqrt(lix*lix + liy*liy);

		// update l0 (either void or contact)
		if (zv[gi] > 0 && zv[ip1[gi]] > 0)
			l0[gi] += dt*(1.0-aL)*cL*(li - l0[gi]);
		else if (zv[gi] <= 0 && NVTOT + growthVerts.size() < NVMAX && li > (1+dt*aL*cL)*l0[gi])
			growthVerts.push_back(gi);

		// -- age angles

		// angle trig functions
		sini = lim1x*liy - lim1y*lix;
		cosi = lim1x*lix + lim1y*liy;

		// angle
		ti = atan2(sini,cosi);

		// update t0
		// if (zv[gi] > 0 || zv[ip1[gi]] > 0 || zv[im1[gi]] > 0)
		// 	t0[gi] = 0.0;
		// else
		t0[gi] += dt*cB*(ti - t0[gi]);


		// age bending mechanical constant
		if (kbi[gi] < kbmax - cKb)
			kbi[gi] += dt*cKb;
	}

	// add vertices
	ngrow = growthVerts.size();
	if (ngrow > 0){
		for (gi=0; gi<growthVerts.size(); gi++){
			cout << "adding vertex between gi=" << growthVerts.at(gi)+gi << " and gip1=" << ip1[growthVerts.at(gi)+gi] << endl;
			addVertex(growthVerts.at(gi)+gi,(1 + (cL/ngrow))*l0[growthVerts.at(gi)+gi]);
		}
	}
	growthVerts.clear();
}


// age shape parameters by "adding" material between contacts
void meso2D::addMesophyllCellMaterial(double dl0){
	// local variables
	bool gtmp = 0;
	int d, gi, gj, gim1, gip1, ctc_im_jn, ctc_ip1m_kp, ci, vi, cin, vin, cip1p, vip1p;
	double dx, dy, dli, dlim1, di, dmin, lix, liy, meanl, testl;

	// vector to pushback vertices where extra material is added
	vector<int> growthVerts;

	// loop over vertices, decide which to grow
	for (gi=0; gi<NVTOT; gi++){
		// info for vertex gi
		cindices(ci,vi,gi);
		meanl = perimeter(ci)/nv[ci];
		testl = 1.001*meanl;
		gip1 = ip1[gi];
		gim1 = im1[gi];

		// check distance to backward neighbor
		dx = x[NDIM*gi] - x[NDIM*gim1];
		dy = x[NDIM*gi + 1] - x[NDIM*gim1 + 1];
		if (pbc[0])
			dx -= L[0]*round(dx/L[0]);
		if (pbc[1])
			dy -= L[1]*round(dy/L[1]);
		dlim1 = sqrt(dx*dx + dy*dy);

		// check distance to forward neighbor
		dx = x[NDIM*gip1] - x[NDIM*gi];
		dy = x[NDIM*gip1 + 1] - x[NDIM*gi + 1];
		if (pbc[0])
			dx -= L[0]*round(dx/L[0]);
		if (pbc[1])
			dy -= L[1]*round(dy/L[1]);
		dli = sqrt(dx*dx + dy*dy);

		// add if neighbors connected to different cells, or if center not connected,
		// but ip1 and im1 are connected
		if (zv[gi] > 0 && zv[gip1] > 0){
			// check min contact to gi
			dmin = 1e6;
			for (gj=0; gj<NVTOT; gj++){
				if (gj < gi)
					gtmp = gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2];
				else if (gj > gi)
					gtmp = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];

				if (gtmp && gi != gj){
					// distance between j and i in x
					dx = x[NDIM*gj] - x[NDIM*gi];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);

					// distance between j and i in y
					dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// full distance
					di = dx*dx + dy*dy;
					if (di < dmin){
						dmin = di;
						ctc_im_jn = gj;
					}
				}
			}

			// check min contact to gip1
			dmin = 1e6;
			for (gj=0; gj<NVTOT; gj++){
				if (gj < gip1)
					gtmp = gij[NVTOT*gj + gip1 - (gj+1)*(gj+2)/2];
				else if (gj > gip1)
					gtmp = gij[NVTOT*gip1 + gj - (gip1+1)*(gip1+2)/2]; 

				if (gtmp && gip1 != gj){
					// distance between j and i in x
					dx = x[NDIM*gj] - x[NDIM*gip1];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);

					// distance between j and i in y
					dy = x[NDIM*gj + 1] - x[NDIM*gip1 + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// full distance
					di = dx*dx + dy*dy;
					if (di < dmin){
						dmin = di;
						ctc_ip1m_kp = gj;
					}
				}
			}

			// if contacts are on different cells, add vertex between contacts if space available
			cindices(cin,vin,ctc_im_jn);
			cindices(cip1p,vip1p,ctc_ip1m_kp);
			if (cin != cip1p){
				cout << "gi=" << gi << ", cin=" << cin << ", gip1=" << gip1 << ", cip1p=" << cip1p << endl;
				growthVerts.push_back(gi);
			}
		}

		// add if lone vertex between two connected vertices, connected to different cells
		else if (zv[gim1] > 0 && zv[gip1] > 0 && zv[gi] == 0){
			// check min contact to gim1
			dmin = 1e6;
			for (gj=0; gj<NVTOT; gj++){
				if (gj < gim1)
					gtmp = gij[NVTOT*gj + gim1 - (gj+1)*(gj+2)/2];
				else if (gj > gim1)
					gtmp = gij[NVTOT*gim1 + gj - (gim1+1)*(gim1+2)/2];

				if (gtmp && gim1 != gj){
					// distance between j and i in x
					dx = x[NDIM*gj] - x[NDIM*gim1];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);

					// distance between j and i in y
					dy = x[NDIM*gj + 1] - x[NDIM*gim1 + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// full distance
					di = dx*dx + dy*dy;
					if (di < dmin){
						dmin = di;
						ctc_im_jn = gj;
					}
				}
			}

			// check min contact to gip1
			dmin = 1e6;
			for (gj=0; gj<NVTOT; gj++){
				if (gj < gip1)
					gtmp = gij[NVTOT*gj + gip1 - (gj+1)*(gj+2)/2];
				else if (gj > gip1)
					gtmp = gij[NVTOT*gip1 + gj - (gip1+1)*(gip1+2)/2]; 

				if (gtmp && gip1 != gj){
					// distance between j and i in x
					dx = x[NDIM*gj] - x[NDIM*gip1];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);

					// distance between j and i in y
					dy = x[NDIM*gj + 1] - x[NDIM*gip1 + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);

					// full distance
					di = dx*dx + dy*dy;
					if (di < dmin){
						dmin = di;
						ctc_ip1m_kp = gj;
					}
				}
			}

			// if contacts are on different cells, add vertex between contacts
			cindices(cin,vin,ctc_im_jn);
			cindices(cip1p,vip1p,ctc_ip1m_kp);
			if (cin != cip1p){
				growthVerts.push_back(gim1);
				growthVerts.push_back(gi);
			}
		}

		// add if void vertex and over extended
		if (zv[gi] == 0){
			// if distance between two vertices is more than twice the vertex radius, add
			if (dlim1 > 4.0*r[gi])
				growthVerts.push_back(gim1);
			if (dli > 4.0*r[gi])
				growthVerts.push_back(gi);
		}

		// // add vertices between neighbors of new vertices
		// else if (zv[gi] == -1){
		// 	// // relax new vertex toward mean segment length of cell
		// 	l0[gim1] += dl0*(meanl - l0[gim1]);
		// 	l0[gi] += dl0*(meanl - l0[gi]);

		// 	// birth if bwd neighbor is connected
		// 	if (zv[gim1] > 0 && dlim1 > testl)
		// 		growthVerts.push_back(gim1);

		// 	// birth vertex if fwd neighbor is connected
		// 	if (zv[gip1] > 0 && dli > testl)
		// 		growthVerts.push_back(gi);
		// }
	}

	// Loop over growth locations, run growth protocol
	if (growthVerts.size() > 0){
		for (gi=0; gi<growthVerts.size(); gi++){
			cout << "adding vertex between gi=" << growthVerts.at(gi)+gi << " and gip1=" << ip1[growthVerts.at(gi)+gi] << endl;
			addVertex(growthVerts.at(gi)+gi,(1+(dl0/growthVerts.size()))*l0.at(growthVerts.at(gi)+gi));
			// addVertex(growthVerts.at(gi)+gi,l0.at(growthVerts.at(gi)+gi));
		}
	}
	growthVerts.clear();
}


// compute a given number of bonded contacts on a single cell
int meso2D::mesoBondedCTCS(int ci){
	// local variables
	bool bondfnd, isConnected;
	int vi, gi, cj, vj, gj, zctmp;

	// loop over pairs of vertices on pairs of cells
	zctmp = 0;
	for (cj=0; cj<NCELLS; cj++){
		// skip if same as ci
		if (cj == ci)
			continue;

		// only check cell pair if somehow in contact
		if (cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2] > 0){
			// since new cell-cell pair, reset bond found to 0
			bondfnd = 0;

			// loop over vertices on cell ci
			gi = szList[ci];
			for (vi=0; vi<nv[ci]; vi++){

				// loop over vertices on cell cj
				gj = szList[cj];
				for (vj=0; vj<nv[cj]; vj++){
					// boolean for bonded contact
					if (gi < gj)
						isConnected = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];
					else
						isConnected = gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2];

					// check if gi and gj bonded, if so, add to cell-cell bonded count, go to next cell
					if (isConnected && !bondfnd){
						zctmp++;
						bondfnd = 1;
						break;
					}

					// increment vertex label on cell cj
					gj++;
				}
				// if bond found, break, increment cj
				if (bondfnd)
					break;

				// increment vertex label on cell ci
				gi++;
			}
		}
	}
	
	// return number of contacts
	return zctmp;
}


// compute the number of bonded contacts between two cells
int meso2D::mesoBondedPAIRS(int ci, int cj){
	// local variables
	bool isConnected;
	int vi, gi, vj, gj, npairs;

	// loop over pairs of vertices on pairs of cells
	npairs = 0;
	gi = szList[ci];
	for (vi=0; vi<nv[ci]; vi++){
		// loop over vertices on cell cj
		gj = szList[cj];
		for (vj=0; vj<nv[cj]; vj++){
			// boolean for bonded contact (gj > gi always here)
			isConnected = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];

			// check if gi and gj bonded, if so, add to cell-cell bonded count, go to next cell
			if (isConnected)
				npairs++;

			// increment vertex label on cell cj
			gj++;
		}

		// increment vertex label on cell ci
		gi++;
	}
	
	// return number of contacts
	return npairs;
}



// add vertex and update all relevant variables
// NOTE: add that only material can be added between contacts of different cells
// * add man file to start from jammed configuration, saves overhead
void meso2D::addVertex(int gi, double newl0){
	// local variables
	int NVVCTS;
	int d, ci, vi, vj, vip1, vim1, cj, gj, gip1 = ip1[gi];
	int gk, gjold, gkold;
	double dx, dy, dr;

	// temporary gij for merging
	NVVCTS = 0.5*NVTOT*(NVTOT-1);
	vector<bool> gijtmp(NVVCTS,0);
	gijtmp = gij;

	// find cell list items
	cindices(ci,vi,gi);

	// get new vertex location
	dx = x[NDIM*gip1] - x[NDIM*gi];
	dy = x[NDIM*gip1 + 1] - x[NDIM*gi + 1];
	if (pbc[0])
		dx -= L[0]*round(dx/L[0]);
	if (pbc[1])
		dy -= L[1]*round(dy/L[1]);
	dr = sqrt(dx*dx + dy*dy);

	// add memory

	// -- dpm shape parameters
	l0.push_back(0.0);
	t0.push_back(0.0);
	r.push_back(0.0);
	im1.push_back(0);
	ip1.push_back(0);

	// -- dynamical variables
	for (d=0; d<NDIM; d++){
		x.push_back(0.0);
		v.push_back(0.0);
		F.push_back(0.0);
	}

	// linked list parameters
	list.push_back(0);

	// -- parameters specific to mesophyll
	kbi.push_back(0.0);
	zv.push_back(0);	

	// print old info to console
	cout << "old:" << endl;
	cout << "vertex " << gi << "; (" << x[NDIM*gi] << ", " << x[NDIM*gi + 1] << ")" << endl;
	cout << "vertex " << gip1 << "; (" << x[NDIM*gip1] << ", " << x[NDIM*gip1 + 1] << ")" << endl;
	cout << "vertex " << ip1[gip1] << "; (" << x[NDIM*ip1[gip1]] << ", " << x[NDIM*ip1[gip1] + 1] << ")" << endl;
	cout << endl;

	// re-sort vertex-based parameters backwards
	for (gj=NVTOT-1; gj>=gi+1; gj--){
		l0.at(gj+1) = l0.at(gj);
		t0.at(gj+1) = t0.at(gj);
		r.at(gj+1) = r.at(gj);
		for (d=0; d<NDIM; d++){
			x.at(NDIM*(gj+1)+d) = x.at(NDIM*gj + d);
			v.at(NDIM*(gj+1)+d) = v.at(NDIM*gj + d);
			F.at(NDIM*(gj+1)+d) = F.at(NDIM*gj + d);
		}
		kbi.at(gj+1) = kbi.at(gj);
		zv.at(gj+1) = zv.at(gj);
	}

	// add to counters
	NVTOT++;
	vertDOF += NDIM;
	nv.at(ci)++;

	// recompute szList, neighbor lists
	for (cj=1; cj<NCELLS; cj++)
		szList.at(cj) = szList.at(cj-1) + nv.at(cj-1);

	// resort ip1 and im1 based on new number of vertices
	for (cj=0; cj<NCELLS; cj++) {
		// vertex indexing
		for (vj=0; vj<nv.at(cj); vj++) {
			// wrap local indices
			vim1 = (vj - 1 + nv.at(cj)) % nv.at(cj);
			vip1 = (vj + 1) % nv.at(cj);

			// get global wrapped indices
			gj = gindex(cj, vj);
			im1.at(gj) = gindex(cj, vim1);
			ip1.at(gj) = gindex(cj, vip1);
		}
	}

	// resize
	// NEED TO MAP OLD CONTACTS TO NEW CONTACTS
	NVVCTS = 0.5*NVTOT*(NVTOT-1);
	gij.clear();
	gij.resize(NVVCTS);
	fill(gij.begin(), gij.end(), 0);
	for (gj=0; gj<NVTOT; gj++){
		if (gj < gi+1)
			gjold = gj;
		else if (gj > gi+1)
			gjold = gj-1;
		else
			continue;

		for (gk=(gj+1); gk<NVTOT; gk++){
			if (gk < gi+1)
				gkold = gk;
			else if (gk > gi+1)
				gkold = gk-1;
			else
				continue;

			gij.at(NVTOT*gj + gk - (gj+1)*(gj+2)/2) = gijtmp.at((NVTOT-1)*gjold + gkold - (gjold+1)*(gjold+2)/2);
		}
	}

	// add information
	x[NDIM*(gi+1)] = x[NDIM*gi] + 0.5*dx;
	x[NDIM*(gi+1) + 1] = x[NDIM*gi + 1] + 0.5*dy;
	v[NDIM*(gi+1)] = 0.0;
	v[NDIM*(gi+1) + 1] = 0.0;
	F[NDIM*(gi+1)] = 0.0;
	F[NDIM*(gi+1) + 1] = 0.0;
	l0[gi+1] = 0.5*newl0;
	l0[gi] = 0.5*newl0;
	// t0[gi+1] = (2.0*PI)/nv.at(ci);
	t0[gi+1] = t0[gi];
	r[gi+1] = r[gi];
	kbi[gi+1] = kbi[gi];
	zv[gi+1] = 0;

	cout << "new:" << endl;
	cout << "vertex " << gi << "; (" << x[NDIM*gi] << ", " << x[NDIM*gi + 1] << ")" << endl;
	cout << "vertex " << gip1 << "; (" << x[NDIM*gip1] << ", " << x[NDIM*gip1 + 1] << ")" << endl;
	cout << "vertex " << ip1[gip1] << "; (" << x[NDIM*ip1[gip1]] << ", " << x[NDIM*ip1[gip1] + 1] << ")" << endl;
	cout << endl;
}


// set t0 to current angles
void meso2D::t0ToCurrent(){
	// local variables
	int gi;
	double lix, liy, lim1x, lim1y, li, ti, sini, cosi;

	for (gi=0; gi<NVTOT; gi++){	
		// segment from i to ip1
		lix = x[NDIM*ip1[gi]] - x[NDIM*gi];
		if (pbc[0])
			lix -= L[0]*round(lix/L[0]);

		liy = x[NDIM*ip1[gi] + 1] - x[NDIM*gi + 1];
		if (pbc[1])
			liy -= L[1]*round(liy/L[1]);

		// segment from im1 to i
		lim1x = x[NDIM*gi] - x[NDIM*im1[gi]];
		if (pbc[0])
			lim1x -= L[0]*round(lim1x/L[0]);

		lim1y = x[NDIM*gi + 1] - x[NDIM*im1[gi] + 1];
		if (pbc[1])
			lim1y -= L[1]*round(lim1y/L[1]);

		// angle trig functions
		sini = lim1x*liy - lim1y*lix;
		cosi = lim1x*lix + lim1y*liy;

		// angle
		ti = atan2(sini,cosi);

		// update t0
		t0[gi] = ti;
	}
}

void meso2D::t0ToReg(){
	int gi, ci, vi;

	gi=0; 
	for (ci=0; ci<NCELLS; ci++){
		for (vi=0; vi<nv[ci]; vi++){
			t0[gi] = (2.0*PI)/nv[ci];
			gi++;
		}
	}
}

// return meso vertex-vertex contact network
void meso2D::getMesoVVContactNetwork(vector<bool> &gijtmp){
	int ci, cj, gi, gj, vi, vj, ctcidx, NVVCTS;
	double rij, sij, dx, dy;

	// check size'
	NVVCTS = 0.5*NVTOT*(NVTOT-1);
	if (gijtmp.size() != NVVCTS || gijtmp.size() != gij.size()){
		cout << "** ERROR: in getMesoVVContactNetwork, input &gijtmp not initialized. Ending. " << endl;
		exit(1);
	}
	else
		cout << "gijtmp size = " << gijtmp.size() << ", gij size = " << gij.size() << endl;

	// update bonded forces
	for (gi=0; gi<NVTOT; gi++){
		for (gj=gi+1; gj<NVTOT; gj++){
			// contact index
			ctcidx = NVTOT*gi + gj - (gi+1)*(gi+2)/2;

			// only add to connections if non-adjacent
			if (gj != ip1[gi] && gj != im1[gi]){
				// contact distance
				sij = r[gi] + r[gj];

				// get vertex-vertex distance
				dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
				dy -= L[1]*round(dy/L[1]);

				dx = x[NDIM*gj] - x[NDIM*gi];
				dx -= L[0]*round(dx/L[0]);

				// get true distance
				rij = sqrt(dx*dx + dy*dy);

				// compute forces and potential
				if (rij < sij)
					gijtmp[ctcidx] = 1;
				else if (gij[ctcidx])
					gijtmp[ctcidx] = 1;


				// // if bonded, set tmp to 1
				// if (gij[ctcidx]) 
				// 	gijtmp[ctcidx] = 1;
				// else {
				// 	// contact distance
				// 	sij = r[gi] + r[gj];

				// 	// get vertex-vertex distance
				// 	dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
				// 	dy -= L[1]*round(dy/L[1]);

				// 	dx = x[NDIM*gj] - x[NDIM*gi];
				// 	dx -= L[0]*round(dx/L[0]);

				// 	// get true distance
				// 	rij = sqrt(dx*dx + dy*dy);

				// 	// compute forces and potential
				// 	if (rij < sij)
				// 		gijtmp[ctcidx] = 1;
				// }
			}
		}
	}

	// sum total number of contacts
	int totalTmpCtcs = 0, totalBndCtcs = 0, prodSum = 0;
	for (gi=0; gi<NVVCTS; gi++){
		if (gijtmp[gi])
			totalTmpCtcs++;
		if (gij[gi])
			totalBndCtcs++;
		if (gij[gi] && gijtmp[gi])
			prodSum++;

	}
	cout << "** gijtmp has " << totalTmpCtcs << " total contacts, gij has " << totalBndCtcs << ", prodSum = " << prodSum << ", size = " << gijtmp.size() << endl;
}


// compute pressure across periodic boundary
double meso2D::mesoInstantaneousPressure(vector<bool> &gijtmp){
	// local variables
	int i, gi, gj, ci, cj, vi, vj, NVVCTS, nvtmp;
	double fa, fli, flim1, fbi, fbim1, cx, cy, xi, yi;
	double flx, fly, fbx, fby;
	double l0im1, l0i, a0tmp, atmp;
	double dx, dy, da, dli, dlim1, dtim1, dti, dtip1;
	double lim2x, lim2y, lim1x, lim1y, lix, liy, lip1x, lip1y, li, lim1;
	double rim2x, rim2y, rim1x, rim1y, rix, riy, rip1x, rip1y, rip2x, rip2y;
	double nim1x, nim1y, nix, niy, sinim1, sini, sinip1, cosim1, cosi, cosip1;
	double ddtim1, ddti;
	double sij, rij, zij, ftmp, fx, fy;
	double P=0.0, dUdL=0.0;

	// check size
	NVVCTS = 0.5*NVTOT*(NVTOT-1);
	if (gijtmp.size() != NVVCTS || gijtmp.size() != gij.size()){
		cout << "** ERROR: in mesoInstantaneousPressure, input &gijtmp not initialized. Ending. " << endl;
		exit(1);
	}

	// reset forces
	resetForcesAndEnergy();

	// loop over vertices, add to forces + dUdL
	U = 0.0;
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

				// add to U
				U += 0.5*ka*da*da;

				// add to dUdL
				dUdL += ((2.0 * ka * atmp)/(a0tmp * L[0])) * da;

				// shape force parameters
				fa = ka * (da / a0tmp);

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
		flim1 	= kl*(1.0/l0im1);
		fli 	= kl*(1.0/l0i);

		// add to forces
		flx 			= (fli*dli*lix/li) - (flim1*dlim1*lim1x/lim1);
		fly 			= (fli*dli*liy/li) - (flim1*dlim1*lim1y/lim1);
		F[NDIM*gi] 		+= flx;
		F[NDIM*gi + 1] 	+= fly;

		// add to U
		U += 0.5*kl*dli*dli;

		// add to dUdL
		dUdL += ((kl * li)/(L[0] * l0i)) * dli;

		// -- Bending force
		fbi = kbi[gi];
		fbim1 = kbi[im1[gi]];

		if (fbi > 0 || fbim1 > 0){
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
			sinim1 = lim2x*lim1y - lim2y*lim1x;
			cosim1 = lim2x*lim1x + lim2y*lim1y;

			sini = lim1x*liy - lim1y*lix;
			cosi = lim1x*lix + lim1y*liy;

			sinip1 = lix*lip1y - liy*lip1x;
			cosip1 = lix*lip1x + liy*lip1y;

			// get normal vectors
			nim1x = lim1y/lim1;
			nim1y = -lim1x/lim1;

			nix = liy/li;
			niy = -lix/li;

			// get change in angles
			dtim1 = atan2(sinim1,cosim1) - t0[im1[gi]];
			dti = atan2(sini,cosi) - t0[gi];
			dtip1 = atan2(sinip1,cosip1) - t0[ip1[gi]];

			// get delta delta theta's
			ddtim1 = (dti - dtim1)/lim1;
			ddti = (dtip1 - dti)/li;

			// add to force
			fbx 			= fbi*ddti*nix - fbim1*ddtim1*nim1x;
			fby 			= fbi*ddti*niy - fbim1*ddtim1*nim1y;
			F[NDIM*gi] 		+= fbx;
			F[NDIM*gi + 1] 	+= fby;

			// update potential energy
			U += 0.5 * kbi[gi] * (dti * dti);
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




	// compute interaction contribution to forces + dUdL
	for (gi=0; gi<NVTOT; gi++){
		for (gj=gi+1; gj<NVTOT; gj++){
			if (gijtmp[NVTOT*gi + gj - (gi+1)*(gi+2)/2]){

				// get cell indices
				cindices(ci,vi,gi);
				cindices(cj,vj,gj);

				// contact distance
				sij = r[gi] + r[gj];

				// get vertex-vertex distance
				dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
				dy -= L[1]*round(dy/L[1]);

				dx = x[NDIM*gj] - x[NDIM*gi];
				dx -= L[0]*round(dx/L[0]);

				// get true distance
				rij = sqrt(dx*dx + dy*dy);

				// // compute forces and potential
				// if (rij > sij){
				// 	// zij: determines strength of bond attraction
				// 	zij 	= 0.5*(zc[ci] + zc[cj])*ctcdel + 1.0;

				// 	// force scale
				// 	ftmp 	= (kc/zij)*(1 - (rij/sij))*(1.0/sij);
				// 	fx 		= ftmp*(dx/rij);
				// 	fy 		= ftmp*(dy/rij);

				// 	// add to U
				// 	U += 0.5*(kc/zij)*pow(1 - (rij/sij),2.0);
				// 	cout << "** In pressure meas, gi=" << gi << " and gj=" << gj << " are stretched bonds, kc/zij = " << kc / zij << "..." << endl;
				// }
				// else {
					// // forces
					// ftmp 	= kc*(1 - (rij/sij))*(1.0/sij);
					// fx 	 	= ftmp*(dx/rij);
					// fy 		= ftmp*(dy/rij);

					// // add to U
					// U += 0.5*kc*pow(1 - (rij/sij),2.0);
				// }
				// forces
				ftmp 	= kc*(1 - (rij/sij))*(1.0/sij);
				fx 	 	= ftmp*(dx/rij);
				fy 		= ftmp*(dy/rij);

				// add to U
				U += 0.5*kc*pow(1 - (rij/sij),2.0);

				// add to net force dot normal for pressure
				F[NDIM*gi] 		-= fx;
				F[NDIM*gi + 1]	-= fy;

				F[NDIM*gj] 		+= fx;
				F[NDIM*gj + 1] 	+= fy;

				// add to dUdL
				dUdL -= ftmp*(rij/L[0]);
			}
		}
	}


	// compute pressure

	// // internal (virial) contribution
	P = 0.0;
	for (i=0; i<vertDOF; i++)
		P += F[i]*x[i];
	P /= (2.0*L[0]*L[1]);

	// // external contribution
	P -= (0.5*dUdL)/L[0];
	// P = -dUdL/(2.0*L[0]);

	// return P value
	return P;
}





// -- Protocol helpers for void bubbles

// relax centers of small bubbles given given cell positions
void meso2D::relaxSmallBubbles(){
	// local variables
	int bi, ci, vi, gi, i;
	int NTOT, NTCELLS;
	double cxtmp, cytmp, dx, dy, rij, sij, ftmp, fxtmp, fytmp;
	double Ftol = 1e-12;

	// get NTOT and # true cells
	NTOT = a0.size();
	NTCELLS = NTOT - NBUBBLES;

	// arrays for moving void bubbles around
	int bbDof = NDIM*NBUBBLES;
	vector<double> bx(bbDof,0.0);
	vector<double> bv(bbDof,0.0);
	vector<double> bf(bbDof,0.0);
	vector<double> br(NBUBBLES,0.0);
	for (bi=0; bi<NBUBBLES; bi++){
		// get global cell index
		ci = NTCELLS + bi;

		// set degree of freedom to be center of cell
		com2D(ci,cxtmp,cytmp);
		bx.at(NDIM*bi) = cxtmp;
		bx.at(NDIM*bi + 1) = cytmp;

		// set radius
		br.at(bi) = 1.25*sqrt(a0[ci]/PI);
	}

	// arrays for fixed cell obstacles
	vector<double> cx(NDIM*NTCELLS,0.0);
	vector<double> cr(NTCELLS,0.0);
	for (ci=0; ci<NTCELLS; ci++){
		// set degree of freedom to be center of cell
		com2D(ci,cxtmp,cytmp);
		cx.at(NDIM*bi) = cxtmp;
		cx.at(NDIM*bi + 1) = cytmp;

		// set radius
		cr.at(bi) = 1.25*sqrt(a0[ci]/PI);
	}

	// FIRE variables
	double P, fnorm, fcheck, vnorm, alpha, dtmax, dtmin;
	int npPos, npNeg, fireit;
	setdt(0.01);

	// Initialize FIRE variables
	P  			= 0;	
	fnorm 		= 0;
	vnorm 		= 0;
	alpha   	= alpha0;

	dtmax   	= 20.0*dt;
	dtmin   	= 1e-2*dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;


	while ((fcheck > Ftol || fireit < NDELAY) && fireit < itmax){
		// compute P
		P = 0.0;
		for (i=0; i<bbDof; i++)
			P += bv[i]*bf[i];

		// print to console
		if (fireit % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << "		B U B B L E  C E N T E R S 			" << endl;
			cout << " 		F I R E 							" << endl;
			cout << "		M I N I M I Z A T I O N 			" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** fireit 	= " << fireit << endl;
			cout << "	** fcheck 	= " << fcheck << endl;
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
				cerr << "	** ERROR: During initial FIRE minimization, P < 0 for too long, so ending." << endl;
				exit(1);
			}

			// take half step backwards, reset velocities
			for (i=0; i<vertDOF; i++){
				// take half step backwards
				bx[i] -= 0.5*dt*bv[i];

				// reset vertex velocities
				bv[i] = 0.0;
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
		for (i=0; i<bbDof; i++)
			bv[i] += 0.5*dt*bf[i];

		// compute fnorm, vnorm and P
		fnorm = 0.0;
		vnorm = 0.0;
		for (i=0; i<bbDof; i++){
			fnorm 	+= bf[i]*bf[i];
			vnorm 	+= bv[i]*bv[i];
		}
		fnorm = sqrt(fnorm);
		vnorm = sqrt(vnorm);

		// update velocities (s.d. vs inertial dynamics) only if forces are acting
		if (fnorm > 0){
			for (i=0; i<bbDof; i++)
				bv[i] = (1 - alpha)*bv[i] + alpha*(bf[i]/fnorm)*vnorm;
		}

		// VV POSITION UPDATE
		for (i=0; i<bbDof; i++)
			bx[i] += dt*bv[i];

		// update forces (function passed as argument)
		for (bi=0; bi<NBUBBLES; bi++){
			for (ci=0; ci<NTCELLS; ci++){
				// get x distance (assume pbc)
				dx = cx[NDIM*ci] - bx[NDIM*bi];
				dx -= L[0]*round(dx/L[0]);

				dy = cx[NDIM*ci + 1] - bx[NDIM*bi + 1];
				dy -= L[1]*round(dy/L[1]);

				// true distance
				rij = sqrt(dx*dx + dy*dy);

				// contact distance
				sij = br[bi] + cr[ci];

				// check overlap
				if (rij < sij){
					ftmp 				= 0.01*(1 - (rij/sij))/sij;
					fxtmp 				= ftmp*(dx/rij);
					fytmp 				= ftmp*(dy/rij);
					bf[NDIM*bi] 		-= fxtmp;
					bf[NDIM*bi + 1] 	-= fytmp;
				}
			}
		}

		// VV VELOCITY UPDATE #2
		for (i=0; i<bbDof; i++)
			bv[i] += 0.5*bf[i]*dt;

		// update fcheck based on fnorm (= force per degree of freedom)
		fcheck = 0.0;
		for (i=0; i<bbDof; i++)
			fcheck += bf[i]*bf[i];
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
		cout << "Bubble center FIRE minimzation converged!" << endl;
	}
}





/******************************

	M E S O P H Y L L

	H E S S I A N

	C O M P U T A T I O N

*******************************/


// hessian for meso cells with vertex-dependent bending energy
void meso2D::mesoBendingHessian(Eigen::MatrixXd &Hb, Eigen::MatrixXd &Sb){
	// local variables
	int gi, kxm1, kym1, kx, ky, kxp1, kyp1, kxp2, kyp2;
	double lxim1, lyim1, lx, ly, ulx, uly, ltmp, si, ci, th;
	double dtiim1x, dtiim1y, dtiix, dtiiy, dtiip1x, dtiip1y;
	double dtip1ix, dtip1iy, dtip1ip2x, dtip1ip2y, dtip1ip1x, dtip1ip1y;
	double dtim1, dti, dtip1;

	// non-dimensionalization
	double Kb;

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
		si = lxim1*ly - lyim1*lx;
		ci = lxim1*lx + lyim1*ly;
		th = atan2(si,ci);
		dth[gi] = th - t0[gi];
	}

	// loop over vertices
	for (gi=0; gi<NVTOT; gi++){
		// bending energy (can be local)
		Kb = kbi[gi];

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


// hessian for meso cells with fixed contact springs
void meso2D::mesoSpringNetworkHessian(Eigen::MatrixXd &Hbnds, Eigen::MatrixXd &Sbnds){
	// local variables
	int ci, vi, gi, cj, vj, gj, kx, ky, lx, ly;
	double sij, dx, dy, rij, kij, h, uxij, uyij, zij;

	// loop over pairs of connected vertices, compute Hessian
	for (gi=0; gi<NVTOT; gi++){
		kx = NDIM*gi;
		ky = kx + 1;
		for (gj=(gi+1); gj<NVTOT; gj++){
			lx = NDIM*gj;
			ly = lx + 1;
			if(gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2]){
				// contact distance
				sij = r[gi] + r[gj];

				// get cell indices
				cindices(ci,vi,gi);
				cindices(cj,vj,gj);

				// zij: determines strength of bond attraction
				zij = 0.5*(zc[ci] + zc[cj])*ctcdel + 1.0;

				// actual distance
				dx = x[lx] - x[kx];
				dy = x[ly] - x[ky];
				if (pbc[0])
					dx -= L[0]*round(dx/L[0]);
				if (pbc[1])
					dy -= L[1]*round(dy/L[1]);
				rij = sqrt(dx*dx + dy*dy);

				// if bond is extended
				if (rij > sij) {
					// spring constant
					kij = kc / (sij * rij * zij);

					// dimensionless overlap
					h = rij / sij;

					// derivatives of distance w.r.t. coordinates
					uxij = dx / rij;
					uyij = dy / rij;

					// compute stiffness and stress matrices (off diagonal, enforce symmetry in lower triangles)

					// -- stiffness matrix
					Hbnds(kx, lx) = -(kc / (sij * sij * zij)) * (uxij * uxij);
					Hbnds(ky, ly) = -(kc / (sij * sij * zij)) * (uyij * uyij);
					Hbnds(kx, ly) = -(kc / (sij * sij * zij)) * (uxij * uyij);
					Hbnds(ky, lx) = -(kc / (sij * sij * zij)) * (uyij * uxij);

					Hbnds(lx, kx) = Hbnds(kx, lx);
					Hbnds(ly, ky) = Hbnds(ky, ly);
					Hbnds(lx, ky) = Hbnds(ky, lx);
					Hbnds(ly, kx) = Hbnds(kx, ly);

					// -- stress matrix
					Sbnds(kx, lx) = kij * (1.0 - h) * (uyij * uyij);
					Sbnds(ky, ly) = kij * (1.0 - h) * (uxij * uxij);
					Sbnds(kx, ly) = -kij * (1.0 - h) * (uxij * uyij);
					Sbnds(ky, lx) = -kij * (1.0 - h) * (uxij * uyij);

					Sbnds(lx, kx) = Sbnds(kx, lx);
					Sbnds(ly, ky) = Sbnds(ky, ly);
					Sbnds(lx, ky) = Sbnds(ky, lx);
					Sbnds(ly, kx) = Sbnds(kx, ly);

					// add to diagonal, using off diagonals and reciprocity

					// -- stiffness matrix
					Hbnds(kx, kx) -= Hbnds(kx, lx);
					Hbnds(ky, ky) -= Hbnds(ky, ly);
					Hbnds(kx, ky) -= Hbnds(kx, ly);
					Hbnds(ky, kx) -= Hbnds(ky, lx);

					Hbnds(lx, lx) -= Hbnds(kx, lx);
					Hbnds(ly, ly) -= Hbnds(ky, ly);
					Hbnds(lx, ly) -= Hbnds(kx, ly);
					Hbnds(ly, lx) -= Hbnds(ky, lx);

					// -- stress matrix
					Sbnds(kx, kx) -= Sbnds(kx, lx);
					Sbnds(ky, ky) -= Sbnds(ky, ly);
					Sbnds(kx, ky) -= Sbnds(kx, ly);
					Sbnds(ky, kx) -= Sbnds(ky, lx);

					Sbnds(lx, lx) -= Sbnds(kx, lx);
					Sbnds(ly, ly) -= Sbnds(ky, ly);
					Sbnds(lx, ly) -= Sbnds(kx, ly);
					Sbnds(ly, lx) -= Sbnds(ky, lx);
				}
			}
		}
	}
}


// function to compute full dynamical matrix for meso network snapshot, compute G
void meso2D::mesoDynamicalMatrix(Eigen::MatrixXd &M, Eigen::MatrixXd &H, Eigen::MatrixXd &S){
	// local variables
	int k, l;

	// print something to the console
	cout << "** Computing dynamical matrix for meso configuration ..." << endl;

	// initialize all possible matrices
	Eigen::MatrixXd Ha(vertDOF, vertDOF);  
	Eigen::MatrixXd Sa(vertDOF, vertDOF);  
	Eigen::MatrixXd Hl(vertDOF, vertDOF);  
	Eigen::MatrixXd Sl(vertDOF, vertDOF);  
	Eigen::MatrixXd Hb(vertDOF, vertDOF);  
	Eigen::MatrixXd Sb(vertDOF, vertDOF);  
	Eigen::MatrixXd Hvv(vertDOF, vertDOF); 
	Eigen::MatrixXd Svv(vertDOF, vertDOF);
	Eigen::MatrixXd Hbnds(vertDOF, vertDOF); 
	Eigen::MatrixXd Sbnds(vertDOF, vertDOF);

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
			Hbnds(k, l) = 0.0;
			Sbnds(k, l) = 0.0;
		}
	}

	// compute all contributions
	cout << "** ** area DM " << endl;
	dpmAreaHessian2D(Ha,Sa);

	cout << "** ** perimeter DM " << endl;
	dpmPerimeterHessian2D(Hl,Sl);

	cout << "** ** bending DM " << endl;
	mesoBendingHessian(Hb,Sb);

	cout << "** ** repulsive springs DM " << endl;
	dpmRepulsiveHarmonicSprings2D(Hvv,Svv);

	cout << "** ** bonds DM " << endl;
	mesoSpringNetworkHessian(Hbnds,Sbnds);

	// construct full dynamical matrix
	cout << "** ** adding together " << endl;
	for (k = 0; k < vertDOF; k++) {
		for (l = 0; l < vertDOF; l++) {
			H(k,l) = Ha(k,l) + Hl(k,l) + Hb(k,l) + Hvv(k,l) + Hbnds(k,l);
			S(k,l) = -Sa(k,l) - Sl(k,l) - Sb(k,l) - Svv(k,l) - Sbnds(k,l);
			M(k,l) = H(k,l) - S(k,l);
		}
	}
}


// function to use dynamical matrix to compute shear modulus, and prints w/ eigenvalues
void meso2D::mesoPrintLinearResponse(meso2DMemFn forceCall, double Ftol, double P0, double dt0){
	// local variables
	int ci, cj, vi, vj, gi, gj, k, l, kx, ky, lx, ly;
	double zij, sij, dx, dy, rij, uy, pglast, pblast;
	double G=0.0, B=0.0;

	// check that hess object is open
	if (!hessout.is_open()){
		cout << "** ERROR: In mesoPrintLinearResponse, hessout object not yet open. Ending. " << endl;
		exit(1);
	}
	else
		cout << "** In mesoPrintLinearResponse, computing and printing G and B... " << endl;

	// // components of dynamical matrix
	// Eigen::MatrixXd M(vertDOF, vertDOF);  
	// Eigen::MatrixXd H(vertDOF, vertDOF);  
	// Eigen::MatrixXd S(vertDOF, vertDOF);

	// // compute
	// mesoDynamicalMatrix(M,H,S);

	// // compute eigenvalues from matrix, plot
	// Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> meig(M);
	// Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> heig(H);
	// Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> seig(S);

	// // eigenvalues
	// Eigen::VectorXd evals = meig.eigenvalues();

	// // also get shear & bulk modulus using numerical derivatives
	G = numericalShearModulus(forceCall, Ftol, dt0, P0);
	B = numericalBulkModulus(forceCall, Ftol, dt0, P0);

	// print eigenvectors
	hessout << setw(w) << left << "NEWFR" << " " << endl;
	hessout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction2D() << endl;

	// print box sizes
	hessout << setw(w) << left << "BOXSZ";
	hessout << setw(wnum) << setprecision(pnum) << left << L[0];
	hessout << setw(wnum) << setprecision(pnum) << left << L[1];
	hessout << endl;

	// print shear modulus
	hessout << setw(w) << left << "SHRMD";
	hessout << setw(wnum) << setprecision(pnum) << left << G;
	hessout << endl;

	// print shear modulus
	hessout << setw(w) << left << "BLKMD";
	hessout << setw(wnum) << setprecision(pnum) << left << B;
	hessout << endl;

	// // print last stress during measurements
	// hessout << setw(w) << left << "STRSS";
	// hessout << setw(wnum) << setprecision(pnum) << left << pglast;
	// hessout << setw(wnum) << setprecision(pnum) << left << pblast;
	// hessout << endl;

	// // print dynamical matrix eigenvalues
	// hessout << setw(w) << left << "MEVAL";
	// for (k=0; k<vertDOF; k++)
	// 	hessout << evals(k) << " ";
	// hessout << endl;

	// // print stiffness matrix eigenvalues
	// evals = heig.eigenvalues();
	// hessout << setw(w) << left << "HEVAL";
	// for (k=0; k<vertDOF; k++)
	// 	hessout << evals(k) << " ";
	// hessout << endl;

	// // print stiffness matrix eigenvalues
	// evals = seig.eigenvalues();
	// hessout << setw(w) << left << "SEVAL";
	// for (k=0; k<vertDOF; k++)
	// 	hessout << evals(k) << " ";
	// hessout << endl;

	// print end frame
	hessout << setw(w) << left << "ENDFR" << " " << endl;
}

// function to compute shear modulus numerically
double meso2D::numericalShearModulus(meso2DMemFn forceCall, double Ftol, double P0, double dt0){
	// local variables
	int k, NVVCTS;
	double sxyold, sxycurr; 

	// shear strain
	double dgamma = 1e-8;
	double gamma = 0.0;
	int NGAMMA = 5;

	// save shear stress
	vector<double> sxyList(NGAMMA+1,0.0);

	// initialize saved variables
	vector<double> Lsave(NDIM,0.0);
	vector<double> stresssave(NDIM + 1,0.0);
	vector<double> lbsave(NDIM,0.0);
	vector<double> xsave(vertDOF, 0.0);
	vector<double> vsave(vertDOF, 0.0);
	vector<double> Fsave(vertDOF, 0.0);
	vector<double> rsave(vertDOF, 0.0);
	vector<double> l0save(vertDOF, 0.0);
	vector<double> t0save(vertDOF, 0.0);
	vector<double> a0save(vertDOF, 0.0);

	Lsave = L;
	stresssave = stress;
	lbsave = lb;
	xsave = x;
	vsave = v;
	Fsave = F;
	rsave = r;
	l0save = l0;
	t0save = t0;
	a0save = a0;

	// initial stress at 0 strain
	sxyList.at(0) = stress[2];

	// initial vertex-vertex contact network
	NVVCTS = 0.5*NVTOT*(NVTOT-1);
	vector<bool> gijtmp(NVVCTS,0);

	// initialize contact network
	getMesoVVContactNetwork(gijtmp);

	// DEBUG: compare init and second relaxation
	int totBndCtcs_real=0, totBndCtcs_test=0;
	double sxx_real, syy_real, sxy_real, sxx_test, syy_test, sxy_test, Fnrm_real=0.0, Fnrm_test=0.0;
	double Fshape_x_real=0.0, Fshape_x_test=0.0, Fshape_y_real=0.0, Fshape_y_test=0.0;
	double Fctc_x_real=0.0, Fctc_x_test=0.0, Fctc_y_real=0.0, Fctc_y_test=0.0;
	double Fnet_x=0.0, Fnet_y=0.0;

	// // real version
	// for (int ci=0; ci<NCELLS; ci++)
	// 	totBndCtcs_real += mesoBondedCTCS(ci);
	// for (int k=0; k<vertDOF; k++)
	// 	Fnrm_real += F[k]*F[k];
	// Fnrm_real = sqrt(Fnrm_real);

	// cout << "=============================" << endl;
	// cout << "COMPARE FORCES BELOW" << endl;
	// cout << "=============================" << endl << endl << endl;



	// real forces
	resetForcesAndEnergy();
	mesoShapeForces();
	for (int gi=0; gi<NVTOT; gi++){
		Fshape_x_real += F[NDIM*gi];
		Fshape_y_real += F[NDIM*gi + 1];
	}
	mesoNetworkForceUpdate();
	Fnet_x = 0.0; 
	Fnet_y = 0.0;
	for (int gi=0; gi<NVTOT; gi++){
		Fnet_x += F[NDIM*gi];
		Fnet_y += F[NDIM*gi + 1];
	}
	Fctc_x_real = Fnet_x - Fshape_x_real;
	Fctc_y_real = Fnet_y - Fshape_y_real;
	sxx_real = stress[0];
	syy_real = stress[1];
	sxy_real = stress[2];


	// forces from shear strain
	cout << endl << endl << endl;
	resetForcesAndEnergy();
	mesoShapeForces(0.0);
	for (int gi=0; gi<NVTOT; gi++){
		Fshape_x_test += F[NDIM*gi];
		Fshape_y_test += F[NDIM*gi + 1];
	}
	mesoNetworkForceUpdate(0.0,gijtmp);
	Fnet_x = 0.0; 
	Fnet_y = 0.0;
	for (int gi=0; gi<NVTOT; gi++){
		Fnet_x += F[NDIM*gi];
		Fnet_y += F[NDIM*gi + 1];
	}
	Fctc_x_test = Fnet_x - Fshape_x_test;
	Fctc_y_test = Fnet_y - Fshape_y_test;
	sxx_test = stress[0];
	syy_test = stress[1];
	sxy_test = stress[2];


	// mesoShearStrainEnthalpyFIRE(0.0, Ftol, P0, dt0, gijtmp);
	// for (int ci=0; ci<NCELLS; ci++)
	// 	totBndCtcs_test += mesoBondedCTCS(ci);
	// for (int k=0; k<vertDOF; k++)
	// 	Fnrm_test += F[k]*F[k];
	// Fnrm_test = sqrt(Fnrm_test);
	

	cout << " * * sxx_real = " << sxx_real << ",  sxx_test = " << sxx_test << endl;
	cout << " * * syy_real = " << syy_real << ",  syy_test = " << syy_test << endl;
	cout << " * * sxy_real = " << sxy_real << ",  sxy_test = " << sxy_test << endl;
	// cout << " * * totBndCtcs_real = " << totBndCtcs_real << ",  totBndCtcs_test = " << totBndCtcs_test << endl;
	// cout << " * * Fnrm_real = " << Fnrm_real << ",  Fnrm_test = " << Fnrm_test << endl;
	cout << " * * Fshape_x_real = " << Fshape_x_real << ", Fshape_x_test = " << Fshape_x_test << endl;
	cout << " * * Fshape_y_real = " << Fshape_y_real << ", Fshape_y_test = " << Fshape_y_test << endl;
	cout << " * * Fctc_x_real = " << Fctc_x_real << ", Fctc_x_test = " << Fctc_x_test << endl;
	cout << " * * Fctc_y_real = " << Fctc_y_real << ", Fctc_y_test = " << Fctc_y_test << endl;

	// cout << "=============================" << endl;
	// cout << "COMPARE FORCES ABOVE" << endl;
	// cout << "=============================" << endl << endl << endl;


	// loop over shear strains gamma, relax at constant pressure using enthalpy FIRE, compute change in Sxy
	for (k=0; k<NGAMMA; k++){
		// update gamma for this iteration
		gamma += dgamma;

		// relax at fixed shear strain + volume, FIXED CONTACT NETWORK
		mesoShearStrainEnthalpyFIRE(gamma, Ftol, P0, dt0, gijtmp);

		// save shear stress
		sxyList.at(k+1) = stress[2];
	}

	// compute numerical derivative, take average
	double Gtmp;
	double G = 0.0;
	gamma = dgamma;
	for (k=0; k<NGAMMA-1; k++){
		Gtmp = -0.5*(sxyList.at(k+2) - sxyList.at(k))/dgamma;
		G += Gtmp;
		cout << " * * k = " << k << ", gamma = " << gamma << ", Gtmp = " << Gtmp << endl;
		gamma += dgamma;
	}
	G /= (NGAMMA-1);
	cout << "on average, G = " << G << endl;

	// reset system
	L = Lsave;
	stress = stresssave;
	lb = lbsave;
	x = xsave;
	v = vsave;
	F = Fsave;
	r = rsave;
	l0 = l0save;
	t0 = t0save;
	a0 = a0save;

	// return numerical result
	return G;
}

// function to compute bulk modulus numerically
// NOTE: need to edit gamma to be strain to CONSTANT PRESSURE!
double meso2D::numericalBulkModulus(meso2DMemFn forceCall, double Ftol, double P0, double dt0){
	// local variables
	int k, gi, d, NVVCTS;
	double Ptmp;

	// shear strain
	double dgamma = 1e-4;
	double gamma = 0.0;
	int NGAMMA = 5;

	// save shear stress
	vector<double> VList(NGAMMA+1,0.0);

	// initialize saved variables
	vector<double> Lsave(NDIM,0.0);
	vector<double> stresssave(NDIM + 1,0.0);
	vector<double> lbsave(NDIM,0.0);
	vector<double> xsave(vertDOF, 0.0);
	vector<double> vsave(vertDOF, 0.0);
	vector<double> Fsave(vertDOF, 0.0);
	vector<double> rsave(vertDOF, 0.0);
	vector<double> l0save(vertDOF, 0.0);
	vector<double> t0save(vertDOF, 0.0);
	vector<double> a0save(vertDOF, 0.0);

	Lsave = L;
	stresssave = stress;
	lbsave = lb;
	xsave = x;
	vsave = v;
	Fsave = F;
	rsave = r;
	l0save = l0;
	t0save = t0;
	a0save = a0;

	// initial stress values at 0 pressure strain
	VList.at(0) = L[0]*L[1];

	// initial vertex-vertex contact network
	NVVCTS = 0.5*NVTOT*(NVTOT-1);
	vector<bool> gijtmp(NVVCTS,0);

	// initialize contact network
	getMesoVVContactNetwork(gijtmp);

	// loop over shear strains gamma, relax using FIRE, compute change in Sxy
	for (k=0; k<NGAMMA; k++){
		// update gamma for this iteration
		gamma += dgamma;

		// update temporary pressure based on gamma
		Ptmp = P0*(1 + gamma);

		// relax at fixed shear strain + volume, FIXED CONTACT NETWORK
		mesoShearStrainEnthalpyFIRE(0.0, Ftol, Ptmp, dt0, gijtmp);

		// save total potential energy + box 
		VList.at(k+1) = L[0]*L[1];
	}

	// compute numerical derivative, take average
	double B = 0.0, Btmp = 0.0;
	double dP, Vtmp, V1, V2, dV, Vavg;
	for (k=0; k<NGAMMA-1; k++){
		// pressure difference between two points
		dP = P0*dgamma;

		// get volume information
		V1 = VList.at(k);
		V2 = VList.at(k+1);
		dV = V2 - V1;
		Vavg = 0.5*(V1 + V2);

		// approximate B
		Btmp = -Vavg*(dP/dV);
		B += Btmp;

		// print
		cout << " * * k = " << k << ", V = " << V1 << ", Btmp = " << Btmp << endl;
	}
	B /= (NGAMMA-1);
	cout << "on average, B = " << B << endl;

	// reset system
	L = Lsave;
	stress = stresssave;
	lb = lbsave;
	x = xsave;
	v = vsave;
	F = Fsave;
	r = rsave;
	l0 = l0save;
	t0 = t0save;
	a0 = a0save;

	// return numerical result
	return B;
}



/******************************

	M E S O P H Y L L

	P R I N T I N G

	F U N C T I O N S

*******************************/


void meso2D::printMesoNetwork2D(){
	// local variables
	bool gtmp;
	int ci, cj, vi, gi, gj, gijtmp, ctmp, zcc, zctmp, zvtmp, zg;
	double xi, yi, dx, dy, Lx, Ly;

	// check if pos object is open
	if (!posout.is_open()){
		cerr << "** ERROR: in printMesoNetwork2D, posout is not open, but function call will try to use. Ending here." << endl;
		exit(1);
	}
	else
		cout << "** In printMesoNetwork2D, printing particle positions to file..." << endl;

	// save box sizes
	Lx = L.at(0);
	Ly = L.at(1);

	// print information starting information
	posout << setw(w) << left << "NEWFR" << " " << endl;
	posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << endl;
	posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction2D() << endl;

	// // print contact info
	// posout << setw(w) << left << "CTCTS";
	// posout << setw(w) << left;
	// for (gi=0; gi<NVTOT; gi++){
	// 	for (gj=(gi+1); gj<NVTOT; gj++){
	// 		gijtmp = NVTOT*gi + gj - (gi+1)*(gi+2)/2;
	// 		if (gij[gijtmp])
	// 			posout << gijtmp << "  ";
	// 	}
	// }

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
		// compute number of contacts with other cells
		zctmp = 0;
		for (cj=0; cj<NCELLS; cj++){
			// compute # of cell contacts
			ctmp = 0;
			if (ci > cj)
				ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2];
			else if (ci < cj)
				ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2];

			if (ctmp > 0)
				zctmp++;
		}

		// cell information
		posout << setw(w) << left << "CINFO";
		posout << setw(w) << left << nv.at(ci);
		posout << setw(w) << left << zctmp;
		posout << setw(wnum) << left << a0.at(ci);
		posout << setw(wnum) << left << area(ci);
		posout << setw(wnum) << left << perimeter(ci);
		posout << endl;

		// get initial vertex positions
		gi = gindex(ci,0);
		xi = x.at(NDIM*gi);
		yi = x.at(NDIM*gi + 1);

		posout << setw(w) << left << "VINFO";
		posout << setw(w) << left << ci;
		posout << setw(w) << left << 0;

		// output initial vertex information
		posout << setw(wnum) << setprecision(pnum) << right << xi;
		posout << setw(wnum) << setprecision(pnum) << right << yi;
		posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << kbi.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << zv.at(gi);
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
			posout << setw(wnum) << setprecision(pnum) << right << kbi.at(gi);
			posout << setw(wnum) << setprecision(pnum) << right << zv.at(gi);
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR" << endl;
}

void meso2D::printMesoNetworkCTCS2D(){
	// local variables
	bool gtmp;
	int ci, cj, vi, gi, gj, gijtmp, ctmp, zcc, zctmp, zvtmp, zg;
	double xi, yi, dx, dy, Lx, Ly;

	// check if pos object is open
	if (!posout.is_open()){
		cerr << "** ERROR: in printMesoNetworkCTCS2D, posout is not open, but function call will try to use. Ending here." << endl;
		exit(1);
	}
	else
		cout << "** In printMesoNetworkCTCS2D, printing particle positions to file..." << endl;

	// save box sizes
	Lx = L.at(0);
	Ly = L.at(1);

	// print information starting information
	posout << setw(w) << left << "NEWFR" << " " << endl;
	posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << endl;
	posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction2D() << endl;

	// print contact info
	posout << setw(w) << left << "CTCTS";
	for (gi=0; gi<NVTOT; gi++){
		for (gj=(gi+1); gj<NVTOT; gj++){
			gijtmp = NVTOT*gi + gj - (gi+1)*(gi+2)/2;
			if (gij[gijtmp])
				posout << gijtmp << "  ";
		}
	}
	posout << endl;

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

	// print instaneous pressure info
	posout << setw(w) << left << "PRESS";
	posout << setw(wnum) << setprecision(pnum) << left << Pinst;
	posout << endl;

	// print coordinate for rest of the cells
	for (ci=0; ci<NCELLS; ci++){
		// compute number of contacts with other cells
		zctmp = 0;
		for (cj=0; cj<NCELLS; cj++){
			// compute # of cell contacts
			ctmp = 0;
			if (ci > cj)
				ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2];
			else if (ci < cj)
				ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2];

			if (ctmp > 0)
				zctmp++;
		}

		// cell information
		posout << setw(w) << left << "CINFO";
		posout << setw(w) << left << nv.at(ci);
		posout << setw(w) << left << zctmp;
		posout << setw(wnum) << left << a0.at(ci);
		posout << setw(wnum) << left << area(ci);
		posout << setw(wnum) << left << perimeter(ci);
		posout << endl;

		// get initial vertex positions
		gi = gindex(ci,0);
		xi = x.at(NDIM*gi);
		yi = x.at(NDIM*gi + 1);

		posout << setw(w) << left << "VINFO";
		posout << setw(w) << left << ci;
		posout << setw(w) << left << 0;

		// output initial vertex information
		posout << setw(wnum) << setprecision(pnum) << right << xi;
		posout << setw(wnum) << setprecision(pnum) << right << yi;
		posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << kbi.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << zv.at(gi);
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
			posout << setw(wnum) << setprecision(pnum) << right << kbi.at(gi);
			posout << setw(wnum) << setprecision(pnum) << right << zv.at(gi);
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR" << endl;
}


void meso2D::printMesoPin2D(vector<double>& xpin, double h){
	// local variables
	bool gtmp;
	int ci, cj, vi, gi, gj, ctmp, zctmp, zvtmp, zg;
	double xi, yi, dx, dy, Lx, Ly;

	// check if pos object is open
	if (!posout.is_open()){
		cerr << "** ERROR: in printMesoPin2D, posout is not open, but function call will try to use. Ending here." << endl;
		exit(1);
	}
	else
		cout << "** In printMesoPin2D, printing particle positions to file..." << endl;

	// save box sizes
	Lx = L.at(0);
	Ly = L.at(1);

	// print information starting information
	posout << setw(w) << left << "NEWFR" << " " << endl;
	posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << endl;
	posout << setw(w) << left << "HPULL" << setw(wnum) << setprecision(pnum) << left << h << endl;

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
		// cell information
		posout << setw(w) << left << "CINFO";
		posout << setw(wnum) << left << xpin[NDIM*ci];
		posout << setw(wnum) << left << xpin[NDIM*ci + 1];
		posout << setw(w) << left << nv.at(ci);
		posout << setw(w) << left << zc.at(ci);
		posout << setw(wnum) << left << a0.at(ci);
		posout << setw(wnum) << left << area(ci);
		posout << setw(wnum) << left << perimeter(ci);
		posout << endl;


		// get initial vertex positions
		gi = gindex(ci,0);
		xi = x.at(NDIM*gi);
		yi = x.at(NDIM*gi + 1);

		posout << setw(w) << left << "VINFO";
		posout << setw(w) << left << ci;
		posout << setw(w) << left << 0;

		// output initial vertex information
		posout << setw(wnum) << setprecision(pnum) << right << xi;
		posout << setw(wnum) << setprecision(pnum) << right << yi;
		posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << kbi.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << zv.at(gi);
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
			posout << setw(wnum) << setprecision(pnum) << right << kbi.at(gi);
			posout << setw(wnum) << setprecision(pnum) << right << zv.at(gi);
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR" << " " << endl;
}


void meso2D::printMesoBondNetwork(){
	// check if pos object is open
	if (!ctcout.is_open()){
		cerr << "** ERROR: in printMesoBondNetwork, ctcout is not open, but function call will try to use. Ending here." << endl;
		exit(1);
	}
	else
		cout << "** In printMesoBondNetwork, printing bond network to file..." << endl;

	// local variables
	int gi, gj;
	for (gi=0; gi<NVTOT; gi++){
		for (gj=(gi+1); gj<NVTOT; gj++){
			ctcout << gij.at(NVTOT*gi + gj - (gi+1)*(gi+2)/2) << " ";
		}
	}
	ctcout << endl;
}

void meso2D::printMesoShearConfigCTCS2D(double gamma){
	// local variables
	bool gtmp;
	int im, ci, cj, vi, gi, gj, gijtmp, ctmp, zcc, zctmp, zvtmp, zg;
	double xi, yi, dx, dy, Lx, Ly;

	// check if pos object is open
	if (!posout.is_open()){
		cerr << "** ERROR: in printMesoShearConfigCTCS2D, posout is not open, but function call will try to use. Ending here." << endl;
		exit(1);
	}
	else
		cout << "** In printMesoShearConfigCTCS2D, gamma = " << gamma << "; printing particle positions & ctcs to file..." << endl;

	// save box sizes
	Lx = L.at(0);
	Ly = L.at(1);

	// print information starting information
	posout << setw(w) << left << "NEWFR" << " " << endl;
	posout << setw(w) << left << "NUMCL" << setw(w) << left << NCELLS << endl;
	posout << setw(w) << left << "PACKF" << setw(wnum) << setprecision(pnum) << left << vertexPackingFraction2D() << endl;
	posout << setw(w) << left << "GAMMA" << setw(wnum) << setprecision(pnum) << left << gamma << endl;

	// print contact info
	posout << setw(w) << left << "CTCTS";
	for (gi=0; gi<NVTOT; gi++){
		for (gj=(gi+1); gj<NVTOT; gj++){
			gijtmp = NVTOT*gi + gj - (gi+1)*(gi+2)/2;
			if (gij[gijtmp])
				posout << gijtmp << "  ";
		}
	}
	posout << endl;

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
		// compute number of contacts with other cells
		zctmp = 0;
		for (cj=0; cj<NCELLS; cj++){
			// compute # of cell contacts
			ctmp = 0;
			if (ci > cj)
				ctmp = cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2];
			else if (ci < cj)
				ctmp = cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2];

			if (ctmp > 0)
				zctmp++;
		}

		// cell information
		posout << setw(w) << left << "CINFO";
		posout << setw(w) << left << nv.at(ci);
		posout << setw(w) << left << zctmp;
		posout << setw(wnum) << left << a0.at(ci);
		posout << setw(wnum) << left << area(ci);
		posout << setw(wnum) << left << perimeter(ci);
		posout << endl;

		// get initial vertex positions
		gi = gindex(ci,0);
		xi = x.at(NDIM*gi);
		yi = x.at(NDIM*gi + 1);

		posout << setw(w) << left << "VINFO";
		posout << setw(w) << left << ci;
		posout << setw(w) << left << 0;

		// output initial vertex information
		posout << setw(wnum) << setprecision(pnum) << right << xi;
		posout << setw(wnum) << setprecision(pnum) << right << yi;
		posout << setw(wnum) << setprecision(pnum) << right << r.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << l0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << t0.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << kbi.at(gi);
		posout << setw(wnum) << setprecision(pnum) << right << zv.at(gi);
		posout << endl;

		// vertex information for next vertices
		for (vi=1; vi<nv.at(ci); vi++){
			// get global vertex index for next vertex
			gi++;

			dy = x.at(NDIM*gi + 1) - yi;
			im = round(dy/Ly);
			dy -= Ly*im;
			yi += dy;

			// get next vertex positions
			dx = x.at(NDIM*gi) - xi;
			dx -= Ly*im*gamma;
			dx -= Lx*round(dx/Lx);
			xi += dx;

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
			posout << setw(wnum) << setprecision(pnum) << right << kbi.at(gi);
			posout << setw(wnum) << setprecision(pnum) << right << zv.at(gi);
			posout << endl;
		}
	}

	// print end frame
	posout << setw(w) << left << "ENDFR" << endl;
}

