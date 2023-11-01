#include "TissueFluid2D.h"

using namespace std;

// constructor
TissueFluid2D::TissueFluid2D(string &inputFileStr, int seed, double boxLengthScale) : adcm2D(inputFileStr, seed, boxLengthScale) {
    // initialize fluid labels 
    fluid_element_label.resize(NVTOT);
    fill(fluid_element_label.begin(), fluid_element_label.end(), -1);
}








// -- Constraints on single cell areas


// function to update areas based on constrained areas
void TissueFluid2D::updateConstrainedAreas(){
    // local variables
    int ci, vi, gi, d;
    double LM;

    // loop over cells and vertices, project areas onto constraints
    gi = 0;
    for (ci=0; ci<NCELLS; ci++){
        // get Lagrange Multiplier
        LM = areaConstraintLM(ci);

        // update positions in normal direction
        for (vi=0; vi<nv[ci]; vi++){
            // update positions
            for (d=0; d<NDIM; d++)
                x[NDIM * gi + d] += LM * normalX(gi, d);

            // update global index
            gi++;
        }
    }
}

// function to update coordinates for cell ci, return lagrange multiplier
double TissueFluid2D::areaConstraintLM(const int ci){
    // local variables
    int vi, gi, nvtmp=nv.at(ci), gip1, gip2;
    double da, C1=0., C2=0., LM;
    double xi, yi, xip1, yip1, xip2, yip2;
    double nxi, nyi, nxip1, nyip1;

    // global index
    gi = szList.at(ci);

    // compute difference in area 
    da = area(ci) - a0.at(ci);

    // precompute first normal
    nxi = normalX(gi, 0);
    nyi = normalX(gi, 1);

    // loop over vectors
    for (vi=0; vi<nvtmp; vi++){
        // relevant indices
        gip1 = ip1[gi];
        gip2 = ip1[gip1];

        xi = relX(gi, 0);
        yi = relX(gi, 1);

        xip1 = relX(gip1, 0);
        yip1 = relX(gip1, 1);

        xip2 = relX(gip2, 0);
        yip2 = relX(gip2, 1);

        // compute normal components
        nxip1 = normalX(gip1, 0);
        nyip1 = normalX(gip1, 1);

        // C1 (linear in normal vector), C2 (quadratic in normal vector)
        C1 += 0.5 * (xi * nyip1 - xip1 * nyi + yip1 * nxi - yi * nxip1);
        C2 += 0.5 * (nxi * nyip1 - nyi * nxip1);
        
        // increment global vertex index
        gi++;

        // update normal for next time
        nxi = nxip1;
        nyi = nyip1;
    }

    // compute lagrange multiplier
    LM = 0.5 * (sqrt((C1 * C1) - (4. * C2 * da)) - C1) / C2;

    // return value
    return LM;
}

// compute forces due to surface tension
void TissueFluid2D::surfaceTensionForceUpdate(){
    // local variables
    int gi, ci, vi;
    double uim1x, uim1y, uix, uiy, li;
    double sti, stim1;

    // loop over all vertices, compute forces due to surface tension
    gi = 0;
    for (ci=0; ci<NCELLS; ci++){
        // initialize segment unit vectors for cell ci
        uim1x = unitSegX(im1[gi], 0);
        uim1y = unitSegX(im1[gi], 1);
        stim1 = max(st[im1[gi]], 0.0);

        // loop over vertices on cell ci
        for (vi=0; vi<nv[ci]; vi++){
            // compute unit vector, surface tension at gi
            li = seg(gi);
            uix = segX(gi, 0) / li;
            uiy = segX(gi, 1) / li;
            sti = max(st[gi], 0.0);

            // update forces
            F[NDIM * gi]        += sti * uix - stim1 * uim1x;
            F[NDIM * gi + 1]    += sti * uiy - stim1 * uim1y;

            // update surface tension, segment geometry
            uim1x = uix;
            uim1y = uiy;
            stim1 = sti;

            // add to potential energy
            U += st[gi] * li;

            // add to pressure and shear stress
            // NOTE: this is not yet the pressure or shear stress, 
            // just the contribution from \partial U / \partial L
            Pinst += (sti * li) / L[0];
            Sinst -= sti * uix * uiy * li;

            // update gi
            gi++;
        }
    }
}