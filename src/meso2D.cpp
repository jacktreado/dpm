/*

	FUNCTION DEFINITIONS for MESO class

	Jack Treado, 06/09/21

*/

#include "meso2D.h"

// namespace
using namespace Eigen;
using namespace std;




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
}


void meso2D::initializeVertexContactNetwork(){
	// local variables
	int NVCTS, gi;

	// check that vertDOF has been assigned
	if (NVTOT <= 0){
		cout << "	** ERROR: in meso2D::initializeVertexContactNetwork, NVTOT not assigned. Ending here." << endl;
		exit(1);
	}
	if (vertDOF <= 0){
		cout << "	** ERROR: in meso2D::initializeVertexContactNetwork, vertDOF not assigned. Ending here." << endl;
		exit(1);
	}
	else if (nv.size() == 0){
		cout << "	** ERROR: in meso2D::initializeVertexContactNetwork, nv vector not assigned. Ending here." << endl;
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
	fill(z.begin(), z.end(), 0);

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
							// get cell and vert indices
							cindices(ci,vi,gi);
							cindices(cj,vj,gj);

							// add bond to network
							if (gi > gj)
								gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
							else if (gi < gj)
								gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1; 

							// update contacts per cell
							z[ci]++;
							z[cj]++;
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
								// get cell and vert indices
								cindices(ci,vi,gi);
								cindices(cj,vj,gj);

								// add bond
								if (gi > gj)
									gij[NVTOT*gj + gi - (gj+1)*(gj+2)/2] = 1;
								else if (gi < gj)
									gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1; 

								// update contacts per cell
								z[ci]++;
								z[cj]++;
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
void meso2D::mesoShapeForces2D(){
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
		fb = kbi[gi]/rho0;

		if (fb > 0){
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


// update forces for bonded mesophyll cells
// * Shape forces of individual cells
// * Repulsive interactions between overlapping vertices
// * Interactions between bonded vertices
void meso2D::mesoNetworkForceUpdate(){
	// local variables
	int gi, gj, ci, cj, vi, vj;
	double rij, sij, zij, dx, dy, fx, fy, ftmp, rho0;

	// normal update (shape + repulsive forces) from base class
	resetForcesAndEnergy();
	mesoShapeForces2D();
	vertexRepulsiveForces2D();

	// update bonded forces
	rho0 = sqrt(a0[0]);
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
					zij = 0.5*(z[ci] + z[cj]) + 1.0;

					// force scale
					ftmp 				= (kc/zij)*(1 - (rij/sij))*(rho0/sij);
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
					stress[0] += (dx*fx*rho0)/(L[0]*L[1]);
					stress[1] += (dy*fy*rho0)/(L[0]*L[1]);
					stress[2] += (0.5*(dx*fy + dy*fx)*rho0)/(L[0]*L[1]);
				}
			}
		}
	}
}


// update forces for mesophyll cells with pin force (ASSUME PIN ALREADY ASSIGNED)
void meso2D::mesoPinForceUpdate(vector<double>& xpin, double kcspring){
	// local variables
	int gi, ci, vi, nvtmp;
	double cx, cy, dcx, dcy;

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

		// add to forces 
		gi = szList[ci];
		for (vi=0; vi<nvtmp; vi++){
			F[NDIM*gi] -= (kcspring/nvtmp)*dcx;
			F[NDIM*gi + 1] -= (kcspring/nvtmp)*dcy;
			gi++;
		}
	}
}



/******************************

	M E S O P H Y L L

	I N T E G R A T O R S

*******************************/

// use FIRE to minimize potential energy in mesophyll cell network
void meso2D::mesoNetworkFIRE(double Ftol, double dt0){
	// local variables
	int i;
	double rho0;

	// check to see if cell linked-list has been initialized
	if (NBX == -1){
		cout << "	** ERROR: In meso2D::mesoNetworkFIRE, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
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

		// update forces in mesophyll network 
		mesoNetworkForceUpdate();

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


// use FIRE to minimize potential energy in mesophyll cells sinking toward center of box
void meso2D::mesoPinFIRE(vector<double>& xpin, double Ftol, double dt0, double kcspring){
	// local variables
	int i, ci, vi, gi, k;
	double cxtmp, cytmp;

	// check to see if cell linked-list has been initialized
	if (NBX == -1){
		cout << "	** ERROR: In meso2D::mesoPinFIRE, NBX = -1, so cell linked-list has not yet been initialized. Ending here.\n";
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

	dtmax   	= 10.0*dt;
	dtmin   	= 1e-2*dt;

	npPos      	= 0;
	npNeg      	= 0;

	fireit    	= 0;
	fcheck  	= 10*Ftol;

	// reset forces and velocities
	resetForcesAndEnergy();
	fill(v.begin(), v.end(), 0.0);

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
void meso2D::mesoNetworkNVE(double T, double dt0, int NT, int NPRINTSKIP){
	// local variables
	int t, i;
	double K;

	// set time step magnitude
	setdt(dt0);

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

		// VV FORCE UPDATE
		mesoNetworkForceUpdate();

		// VV VELOCITY UPDATE #2
		for (i=0; i<vertDOF; i++)
			v[i] += 0.5*F[i]*dt;

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
			if (xtraout.is_open())
				xtraout << t << " " << U << " " << K << endl;
			if (posout.is_open())
				printConfiguration2D();
		}
	}
}











/******************************

	M E S O P H Y L L

	P R O T O C O L S

*******************************/


// simulate network extension using decompression
void meso2D::mesoNetworkExtension(double Ftol, double dt0, double delShrink, double dphiPrint, double phiMin){
	// local variables
	int k = 0;
	double lastPrintPhi, scaleFactor = 1.0 - 2.0*delShrink;

	// current preferred packing fraction
	double phi0 = vertexPreferredPackingFraction2D();
	lastPrintPhi = phi0;

	// loop until phi0 < phiMin
	while (phi0 > phiMin && k < itmax){
		// relax current configuration
		mesoNetworkFIRE(Ftol, dt0);

		// break contact network
		updateMesophyllBondNetwork();

		// age particle shapes
		ageMesophyllShapeParameters();

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
		cout << "	* lastPrintPhi 	= " << lastPrintPhi << endl;
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
		if ((lastPrintPhi - phi0) > dphiPrint){
			// print positions
			printConfiguration2D();

			// update last phi when printed, for next time
			lastPrintPhi = phi0;
		}

		// scale particle sizes
		scaleParticleSizes2D(scaleFactor);

		// update new packing fraction
		phi0 = vertexPreferredPackingFraction2D();

		// update iterate
		k++;
	}
}


// update mesophyll network using MC
void meso2D::updateMesophyllBondNetwork(){
	// local variables
	bool isConnected;
	int ci, cj, vi, vj, gi, gj;
	double dx, dy, sij, rij, dU, poff, h=0.5, h2=h*h, rdraw;

	// loop over pairs of vertices, check whether to connect or detach
	for (gi=0; gi<NVTOT; gi++){
		cindices(ci,vi,gi);
		for (gj=gi+1; gj<NVTOT; gj++){
			cindices(cj,vj,gj);
			if (ci != cj){
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

				// get current contact status
				isConnected = gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2];

				// reconnect bond if vertices come into contact
				if (rij <= sij && !isConnected){
					gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 1;
					z[ci]++;
					z[cj]++;
				}
				else if (rij > sij && isConnected){
					// change in energy from bond breaking
					dU = 1.0 - (pow(1 - (rij/sij),2.0)/h2);

					// remove if bond detaching decreases energy
					if (dU < 0){
						// remove vertex-vertex bond
						gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 0;

						// remove single cell-cell contact
						if (ci > cj)
							cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]--;
						else
							cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]--; 
						z[ci]--;
						z[cj]--;
					}
					else{
						// else, remove conditionally
						poff = exp(-betaEff*dU);
						rdraw = drand48();

						// detach
						if (poff > rdraw){
							// detach vv contact
							gij[NVTOT*gi + gj - (gi+1)*(gi+2)/2] = 0;

							// remove single cell-cell contact
							if (ci > cj)
								cij[NCELLS*cj + ci - (cj+1)*(cj+2)/2]--;
							else
								cij[NCELLS*ci + cj - (ci+1)*(ci+2)/2]--; 
							z[ci]--;
							z[cj]--;
						}
					}
				}
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

		// update l0
		l0[gi] += cL*(li - l0[gi]);

		// -- age angles

		// angle trig functions
		sini = lix*lim1y - liy*lim1x;
		cosi = lix*lim1x + liy*lim1y;

		// angle
		ti = atan2(sini,cosi);

		// update t0
		t0[gi] += cB*(ti - t0[gi]);


		// age bending mechanical constant
		kbi[gi] += cKb;
	}
}



// drag cell pins away from center
void meso2D::mesophyllPinExtension(double Ftol, double dt0, double hmax, double dh, double dhprint, double kcspring){
	// local variables
	int k=0, ci;
	double cx, cy, dcx, dcy, rho0, h=0.0, lastPrinth=0.0;
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

	// loop over extension steps
	rho0 = sqrt(a0[0]);
	while (h < hmax && k < itmax){
		// relax current configuration
		mesoPinFIRE(xpin,Ftol,dt0,kcspring);

		// update contact network
		updateMesophyllBondNetwork();

		// age particle shapes
		ageMesophyllShapeParameters();

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
		for (ci=0; ci<NCELLS; ci++)
			cout << "px=" << xpin[NDIM*ci] << ",   py=" << xpin[NDIM*ci+1] << endl;

		// print positions if change in packing fraction is large enough
		if (abs(lastPrinth - h) > dhprint){
			// print positions
			printConfiguration2D();

			// update last print h
			lastPrinth = h;
		}

		// move pins
		for (ci=0; ci<NCELLS; ci++){
			xpin[NDIM*ci] += dh*rho0*cos(th[ci]);
			xpin[NDIM*ci + 1] += dh*rho0*sin(th[ci]);
		}

		// update new h
		h += dh;

		// update iterate
		k++;
	}
}





















