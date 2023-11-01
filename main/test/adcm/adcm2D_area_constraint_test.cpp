// Function to test adcm2D area constraint function
// compile with: g++ --std=c++14 -I src main/test/adcm/adcm2D_area_constraint_test.cpp src/dpm.cpp src/adcm2D_*.cpp -o adcm2D_area_constraint_test

#include "adcm2D.h"

using namespace std;

int main(){
    // print to console
    cout << "** In test: adcm2D_area_constraint_test.cpp, testing area constraint functionality" << endl;

    // local variables
    int NCELLS = 1;
    int NV = 24;
    double randforce = 0.0, RFMAG = 0.1;
    string fstr = "adcm2D_area_constraint_test.pos";

    cout << "=================================================================" << endl;
    cout << "** TEST 1: make a single particle, one step, check change to area" << endl;
    cout << "=================================================================" << endl;

    // instantiate an ADCM2D system with a single particle
    adcm2D sys1(NCELLS, NV, 0, 0.25, 3.0, 1.0, 0.001, 1);

    // open output file
    sys1.openPosObject(fstr);

    // measure initial area
    double astart = sys1.area(0);

    // make initial perturbation
    srand48(time(NULL));
    for (int i=0; i<NV; i++){
        for (int d=0; d<2; d++){
            // random force
            randforce = RFMAG * drand48();

            // perturb
            sys1.addx(randforce, i, d);
        }
    }

    // check new area
    double aend = sys1.area(0);
        
    // make constraint correction
    sys1.updateConstrainedAreas();

    // measure new area, check to see if constraint worked
    double acorrect = sys1.area(0);

    // loop over vertices, check definition of normals
    int im1, ip1;
    double nix, niy, nixReal, niyReal;
    for (int i=0; i<NV; i++){
        nix = sys1.normalX(i, 0);
        niy = sys1.normalX(i, 1);

        ip1 = (i + 1) % NV;
        im1 = (i + NV - 1) % NV;
        nixReal = 0.5 * (sys1.getx(ip1, 1) - sys1.getx(im1, 1));
        niyReal = 0.5 * (sys1.getx(im1, 0) - sys1.getx(ip1, 0));

        cout << "nix = " << nix << ", nixReal = " << nixReal << ";     niy = " << niy << ", niyReal = " << niyReal << ";   diff = " << nix + niy - nixReal - niyReal << endl;
    }

    // print to console
    cout << "** Area update complete, results:" << endl;
    cout << setprecision(12) << endl;
    cout << " -- astart = " << astart << endl;
    cout << " -- aend = " << aend << endl;
    cout << " -- acorrect = " << acorrect << endl;
    cout << " -- a0 = " << sys1.geta0(0) << endl;


    // -- TEST 2: impose some external pressure by increasing surface tension, make sure

    // reset to center of box
    sys1.setCom(0, 0.5 * sys1.getL(0), 0.5 * sys1.getL(1));

    // set surface tension

    
    // area stays constant during steepest descent minimization
    // while (t < tmax && k < kmax){
    //     // update forces
    //     adcm2DSurfaceTensionForceUpdate();


    //     // update iterate
    //     kmax++;
    // }
    // if (k == kmax)
        // cout << "ERROR in adcm2D_area_constraint_test, k -> kmax = " << k << ", ending." << endl;

}