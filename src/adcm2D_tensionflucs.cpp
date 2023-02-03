/*

	FUNCTION DEFINITIONS for ADCM2D class
    for tension fluctuation code

    * Stands for Active Deformable Cell Model 
	* only for use in two dimensions

	Jack Treado, 09/13/22

*/

#include "adcm2D.h"

// namespacess
using namespace std;


// update all forces when using active tension fluctuations
void adcm2D::activeTensionForceUpdate(){
	// update vertex numbers
	// checkVertices();

	// reset energies, stresses, forces
	U = 0.0; 
	Pinst = 0.0; 
	Sinst = 0.0; 
	fill(F.begin(), F.end(), 0.0);

	// use whatever shape force
	(*this.*shpFrc)();

	// reset surface tension to the void values
	int gi = 0;
	for (int ci=0; ci<NCELLS; ci++){
		for (int vi=0; vi<nv.at(ci); vi++){
			// set tension cell-void coupling
			st.at(gi) = stMat.at(0).at(ci+1);
			// st.at(gi) += dt * (stMat.at(0).at(ci+1) - st.at(gi));

			// increment global index
			gi++;
		}
	}

	// use pairwise force update
	circuloLinePWForceUpdate();

	// update stresses
	stressUpdate();
}

void adcm2D::activeTensionForceUpdateVertexPreservation(){
	// reset energies, stresses, forces
	U = 0.0; 
	Pinst = 0.0; 
	Sinst = 0.0; 
	fill(F.begin(), F.end(), 0.0);

	// use whatever shape force
	(*this.*shpFrc)();

	// reset surface tension to the void values
	int gi = 0;
	for (int ci=0; ci<NCELLS; ci++){
		for (int vi=0; vi<nv.at(ci); vi++){
			// set tension cell-void coupling
			st.at(gi) = stMat.at(0).at(ci+1);
			// st.at(gi) += dt * (stMat.at(0).at(ci+1) - st.at(gi));

			// increment global index
			gi++;
		}
	}

	// use pairwise force update
	circuloLinePWForceUpdate();

	// update stresses
	stressUpdate();
}


// short-range (SR) _attractive_ pairwise force between two circulolines WITH SURFACE TENSION FLUCTUATIONS
// 
// Procedure:
// 1. Compute dlijApprox, approx dist (lower bound) between segment centers
// 2. If dlijApprox small enough, 
// 		a. 	compute projection from edges im1, i, ip1 onto vertex j
// 				i. 		compute vertex-vertex force if vertex i and vertex j overlap
// 				ii.		otherwise, compute edge-vertex force depending on projection overlap between edges on mu and vertex (j, nu)
// 		b. 	compute projection from edges jm1, j, jp1 onto vertex i
// 				i. 		compute only edge-vertex force depending on projection overlap between edges on nu and vertex (i, mu)
// 
// TO DO:
//  * How to make surface surface tensions vary properly?
// 	-- what happens if vertex j overlaps with 1 or more edges near i, but vertex i overlaps also with multiple edges near j? 
// 	-- what happens if the number of edges overlapped by vertex j does not equal the number of edges overlapped by vertex i? 
// 		** will this create problems?
//	-- how to make sure surface tension variation is not extreme when edges bind and unbind? How to lag the surface tensions on a segment-basis across multiple functions?
// 
void adcm2D::SRAttractiveActiveTensionPWForce(const int gi, const int gj, bool &vivj, bool &viej, bool &vjei) {
	// local variables
	int ci, vi, cj, vj;
	double ftmp;

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

	// test values, when dr < sij
	int pi, pj;
	double dij, dijx, dijy, dji, djix, djiy, ti, tj;

	// get ci and cj (for surface tension update)
	cindices(ci, vi, gi);
	cindices(cj, vj, gj);

	// contact distances
	double sij = r[gi] + r[gj];
	double shellij = (1.0 + l2)*sij;
	double cutij = (1.0 + l1)*sij; 

	// indices
	int gip1 = ip1[gi];
	int gim1 = im1[gi];
	int gjp1 = ip1[gj];
	int gjm1 = im1[gj];

	// distance between vertex i and j
	double dx = deltaX(gi, gj, 0);
	double dy = deltaX(gi, gj, 1);
	double dr = sqrt(dx * dx + dy * dy);

	// distance between vertex i and jp1
	double drjip1 = deltaR(gj, gip1);
	double drijp1 = deltaR(gi, gjp1);

	// Gut check distance
	double lij = 0.5 * (seg(gi) + seg(gj));
	double dlijApprox = 0.5 * (dr + deltaR(gip1, gjp1));
	if (dlijApprox < 6.0 * lij) {
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
				if (ti < 0)
					ftmp = vvSoftAdhesionForce(pi, gj, dij, dijx, dijy);
				else
					ftmp = evSoftAdhesionForce(pi, gj, dij, dijx, dijy, ti);
			}
			else {
				if (tj < 0)
					ftmp = vvSoftAdhesionForce(pj, gi, dji, djix, djiy);
				else
					ftmp = evSoftAdhesionForce(pj, gi, dji, djix, djiy, tj);
			}
		}
		else if (drjip1 > shellij && hr_ij < shellij && tij > 0 && tij < 1) {
			// compute projection onto adjacent edges
			hr_im1j = edge2VertexDistance(gim1, gj, hx_im1j, hy_im1j, tim1j);
			hr_ip1j = edge2VertexDistance(gip1, gj, hx_ip1j, hy_ip1j, tip1j);

			// check location and projections to determine force
			if ((tim1j < 0 || tim1j > 1) && (tip1j < 0 || tip1j > 1)){
				// IF vertex gj projects onto "bulk" of edge gi, no other edge
				ftmp = evSoftAdhesionForce(gi, gj, hr_ij, hx_ij, hy_ij, tij);
			}				
			else if (tim1j > 0 && tim1j < 1) {
				if (hr_ij < hr_im1j)
					ftmp = evSoftAdhesionForce(gi, gj, hr_ij, hx_ij, hy_ij, tij);				// IF vertex gj projects onto edges gi and gim1, closer to gi (REGION III')
				else
					ftmp = evSoftAdhesionForce(gim1, gj, hr_im1j, hx_im1j, hy_im1j, tim1j);		// IF vertex gj projects onto edges gi and gim1, closer to gim1 (REGION II')
			}
			else if (tip1j > 0 && tip1j < 1) {
				if (hr_ij < hr_ip1j)
					ftmp = evSoftAdhesionForce(gi, gj, hr_ij, hx_ij, hy_ij, tij);				// IF vertex gj projects onto edges gi and gip1, closer to gi (REGION II' centered on ip1)
				else 
					ftmp = evSoftAdhesionForce(gip1, gj, hr_ip1j, hx_ip1j, hy_ip1j, tip1j);     // IF vertex gj projects onto edges gi and gip1, closer to gip1 (REGION III' centered on ip1)
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
				ftmp = evSoftAdhesionForce(gj, gi, hr_ji, hx_ji, hy_ji, tji);
			}
			else if (tjm1i > 0 && tjm1i < 1){
				if (hr_ji < hr_jm1i)
					ftmp = evSoftAdhesionForce(gj, gi, hr_ji, hx_ji, hy_ji, tji);				// IF vertex gi projects onto edges gj and gjm1, closer to gj (REGION III')
				else
					ftmp = evSoftAdhesionForce(gjm1, gi, hr_jm1i, hx_jm1i, hy_jm1i, tjm1i);		// IF vertex gi projects onto edges gj and gjm1, closer to gjm1 (REGION II')
			}
			else if (tjp1i > 0 && tjp1i < 1) {
				if (hr_ji < hr_jp1i)
					ftmp = evSoftAdhesionForce(gj, gi, hr_ji, hx_ji, hy_ji, tji);				// IF vertex gi projects onto edges gj and gjm1, closer to gj (REGION II' centered on ip1)
				else 
					ftmp = evSoftAdhesionForce(gjp1, gi, hr_jp1i, hx_jp1i, hy_jp1i, tjp1i);		// IF vertex gi projects onto edges gj and gjm1, closer to gjm1 centered on ip1)
			}
		}
	}
}
