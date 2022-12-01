// Main file to measure friction coefficient between two contacting cells
// compile with: g++ -O3 --std=c++11 -I src main/adcm2D/measureFriction.cpp src/*.cpp -o fric.o
// run with: ./fric.o 24 1e-3 0.1 1 pos.test

#include "adcm2D.h"
#include <sstream>

// constants in sim
# define NCELLS 2

using namespace std;

int main(int argc, char const *argv[]){
    // simulation details that are fixed
    const double dt = 0.01;
    const double kc = 1.0;
    const double phi0 = 0.2;
    const double boxLengthScale = 3.0;
    const double clScale = 1.0;
    const double Ftol = 5e-7;
    const int kmax = 1e7;

    // local variables
    int numverts, seed;
    double kctr, gam0;

    // inputs
    string numverts_str = argv[1];
    string gam0_str = argv[2];
    string kctr_str = argv[3];
    string seed_str = argv[4];
    string pos_str = argv[5];

    // convert to sstreams
    stringstream numverts_ss(numverts_str);
    stringstream gam0_ss(gam0_str);
    stringstream kctr_ss(kctr_str);
    stringstream seed_ss(seed_str);

    // convert to numbers
    numverts_ss >> numverts;
    gam0_ss >> gam0;
    kctr_ss >> kctr;
    seed_ss >> seed;

    // instantiate object
    adcm2D sim(NCELLS, numverts, 0.0, phi0, boxLengthScale, clScale, gam0, seed);

    // open output file, print
	sim.openPosObject(pos_str);

    // set simulation constants
    sim.setkc(kc);
    sim.setdt_pure(dt);
    sim.setW(0);
    sim.setl2(1.0);
    sim.setl1(0.1);
    sim.useAttractiveForce();
    // sim.useAlignmentForce();

    // set particles in middle of box
    double L = sim.getL(0);
    sim.setCom(0, 0.25 * L, 0.5 * L);
    sim.setCom(1, 0.75 * L, 0.5 * L);

    // relax cells together using central potential
    

    // -- Protocol to relax cells toward each other
    double fcheck = 10 * Ftol;
    int k = 0;
    int NSKIP = round(100 * sqrt(numverts / kctr));
    // int NSKIP = 50;
    bool touched = false;

    // variables for central potential force computation
    double cx, cy, dcx, dcy, dfx, dfy;
    
    // loop over gradient descent dynamics
    while (k < kmax && fcheck > Ftol){
        // update usual forces
        sim.adcm2DForceUpdate();

        // add central potential
        for (int ci=0; ci<NCELLS; ci++){
            // get centers of mass
            sim.com2D(ci, cx, cy);

            // compute force
            dcx = 0.5 * L - cx;
            dcy = 0.5 * L - cy;

            // add force to center of mass if far enough away
            if (abs(dcx) > 0.4 && !touched){
                dfx = kctr * dcx;
                dfy = kctr * dcy;
                
                sim.addComF(ci, 0, dfx);
                sim.addComF(ci, 1, dfy);                
            }
            else if (abs(dcx) < 0.4) {
                touched = true;
            }
        }

        // measure force RMSD
        fcheck = 0.0;
        for (int gi=0; gi<sim.getNVTOT(); gi++)
            fcheck += pow(sim.getF(gi, 0), 2.0) + pow(sim.getF(gi, 1), 2.0);
        fcheck = sqrt(fcheck / sim.getNVTOT());

        // update positions
        for (int gi=0; gi<sim.getNVTOT(); gi++){
            for (int d=0; d<sim.getNDIM(); d++){
                sim.addx(sim.getdt() * sim.getF(gi, d), gi, d);
            }
        }
            

        // print status
        if (k % NSKIP == 0){
            cout << endl << endl;
            cout << "=== Relaxing toward center === " << endl;
            cout << "* k = " << k << endl;
            cout << "* fcheck = " << fcheck << " / " << Ftol << endl;
            cout << "* dcx = " << dcx << endl;
            cout << "* dcy = " << dcy << endl;
            sim.printADCM2DConfiguration();
        }

        // update iterate
        k++;
    }
}