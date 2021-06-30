/*

	FUNCTION DEFINITIONS for tumor class

	Jack Treado, 06/04/21

	** TO DO 06/24/21
	1. Add adipocyte / tumor boundary sim initialization function
	2. Make protocol that compresses ONLY so that we have premade initial conditions, don't have to spend overhead
	3. Add constant pressure via piston attached to left-side wall, box vol. changes as cells invade
		* Need to add kinetic term to stress tensor

*/


#include "tumor.h"
#include "dpm.h"

// namespace
using namespace Eigen;
using namespace std;



/*********************************

	C O N S T R U C T O R S 

	&

	D E S T R U C T O R S 

**********************************/




/*********************************

	T U M O R   C E L L

	I N I T I A L I Z A T I O N

**********************************/


void tumor::initializeSingleTumorCell(){
	// local variables
	int vi;
	int dtmp;

	// only proceed if ntmp has not been set yet, but nv etc has
	if (ntmp != 0){
		cout << "\t ** ERROR: in initializeSingleTumorCell, ntmp = " << ntmp << " which is not 0, cannot reset. Ending here." << endl;
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

	// initialize ntmp to 1
	ntmp = 1;

	// create first cell centered at origin
	for (vi=0; vi<nv.at(0); vi++){
		// length from center to vertex
		dtmp = sqrt((2.0*a0.at(0))/(nv.at(0)*sin((2.0*PI)/nv.at(0))));

		// set positions
		x.at(NDIM*vi) 		= dtmp*cos((2.0*PI*vi)/nv.at(0)) + 1e-2*l0[vi]*drand48();
		x.at(NDIM*vi + 1)	= dtmp*sin((2.0*PI*vi)/nv.at(0)) + 1e-2*l0[vi]*drand48();
	}
}








/*********************************

	E D I T I N G   &

			U P D A T I N G

**********************************/

// divide a single cell, assume preallocation
void tumor::divide(int ci){
	// local variable
	int gi, g1tmp, g2tmp, vi, ci1, ci2, nv1, nv2, icut1, icut2, nh1, nh2, vitmp;
	double Dx, Dy, Dnorm, dhatx, dhaty, nhatx, nhaty, xtmp, ytmp;

	// cell indices
	ci1 = ntmp - 1;
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

	F O R C E 

	U P D A T E S

*******************************/

void tumor::tumorInteractionForces2D(){
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

				// attraction distances
				shellij = (1.0 + l2)*sij;
				cutij = (1.0 + l1)*sij;

				// particle distance
				dx = x[NDIM*gj] - x[NDIM*gi];
				if (pbc[0])
					dx -= L[0]*round(dx/L[0]);
				if (dx < shellij){
					dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
					if (pbc[1])
						dy -= L[1]*round(dy/L[1]);
					if (dy < shellij){
						rij = sqrt(dx*dx + dy*dy);
						if (rij < shellij){
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
							stress[1] 			+= dx*fx;
							stress[2] 			+= dy*fy;
							stress[3] 			+= 0.5*(dx*fy + dy*fx);

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

					// attraction distances
					shellij = (1.0 + l2)*sij;
					cutij = (1.0 + l1)*sij;

					// particle distance
					dx = x[NDIM*gj] - x[NDIM*gi];
					if (pbc[0])
						dx -= L[0]*round(dx/L[0]);
					if (dx < shellij){
						dy = x[NDIM*gj + 1] - x[NDIM*gi + 1];
						if (pbc[1])
							dy -= L[1]*round(dy/L[1]);
						if (dy < shellij){
							rij = sqrt(dx*dx + dy*dy);
							if (rij < shellij){
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
								stress[1] 			+= dx*fx;
								stress[2] 			+= dy*fy;
								stress[3] 			+= 0.5*(dx*fy + dy*fx);

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








