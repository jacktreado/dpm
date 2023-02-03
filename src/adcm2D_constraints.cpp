/*

	FUNCTION DEFINITIONS for ADCM2D class
    for code with area constraints

	Jack Treado, 02/03/23

*/

#include "adcm2D.h"

// namespacess
using namespace std;



// function to update areas based on constrained areas
void adcm2D::updateConstrainedAreas(){
    // local variables
    int ci, vi, gi, gip1, gip2, gk;
    double da, C1=0., C2=0., LM;
    double xi, yi, xip1, yip1, xip2, yip2;
    double nxi, nyi, nxip1, nyip1;

    // constraint forces
    double gx, gy;

    // store normals
    vector<double> nx;
    vector<double> ny;

    // loop over cells and vertices, project areas onto constraints
    gi = 0;
    for (ci=0; ci<NCELLS; ci++){
        // compute difference in area 
        da = area(ci) - a0.at(ci);

        // reset vector to be nv long
        nx.resize(nv[ci]);
        ny.resize(nv[ci]);

        // precompute first normal
        nxi = 0.5 * (x[NDIM * ip1[gi] + 1] - x[NDIM * im1[gi] + 1]);
        nyi = 0.5 * (x[NDIM * im1[gi]] - x[NDIM * ip1[gi]]);

        // compute constants for cell ci
        for (vi=0; vi<nv[ci]; vi++){
            // relevant indices
            gip1 = ip1[gi];
            gip2 = ip1[gip1];

            xi = x[NDIM * gi];
            yi = x[NDIM * gi + 1];

            xip1 = x[NDIM * gip1];
            yip1 = x[NDIM * gip1 + 1];

            xip2 = x[NDIM * gip2];
            yip2 = x[NDIM * gip2 + 1];

            // compute / store normal components
            nx.at(vi) = nxi;
            ny.at(vi) = nyi;
            nxip1 = 0.5 * (yip2 - yi);
            nyip1 = 0.5 * (xi - xip2);

            // C1 (quadratic in normal vector), C2 (linear in normal vector)
            C2 += nxi * nyip1 - nyi * nxip1;
            C2 += xi * nyip1 - xip1 * nyi + yip1 * nxi - yi * nxip1;
            
            // increment global vertex index
            gi++;

            // update normal for next time
            nxi = nxip1;
            nyi = nyip1;
        }

        // compute lagrange multiplier
        LM = 0.5 * (sqrt(C1 * C1 - 4. * C2 * da) - C1) / C2;

        // loop again, apply constraint
        gk = szList[ci];
        for (vi=0; vi<nv[vi]; vi++){
            // apply constraint
            x[NDIM * gk]        += LM * nx[vi];
            x[NDIM * gk + 1]    += LM * ny[vi];

            // update index for coordinates
            gk++;
        }

        // reset vectors
        nx.clear();
        ny.clear();
    }
    
}