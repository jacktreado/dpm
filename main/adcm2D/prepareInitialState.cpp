// Main file to prepare initial state at positive repulsive pressure
// compile:     g++ --std=c++11 -O3 -I src/ main/adcm2D/prepareInitialState.cpp src/*.cpp -o prep.o
// test:        ./prep.o 24 24 0.1 1e-4 5e-4 1 pos.test

#include "adcm2D.h"
#include <sstream>

using namespace std;

int main(int argc, char const *argv[]){
    // simulation details that are fixed
    const double dt = 0.01;
    const double kc = 1.0;
    const double phi0 = 0.5;
    const double boxLengthScale = 3.0;
    const double clScale = 1.0;
    const double Ftol = 1e-10;
    const double dPtol = 1e-10;

    // local variables
    int NCELLS, NV, seed;
    double disp, P0, gam0;

    // inputs
    string numcells_str = argv[1];
    string numverts_str = argv[2];
    string disp_str = argv[3];
    string gam0_str = argv[4];
    string P0_str = argv[5];
    string seed_str = argv[6];
    string pos_str = argv[7];

    // convert to sstreams
    stringstream numcells_ss(numcells_str);
    stringstream numverts_ss(numverts_str);
    stringstream disp_ss(disp_str);
    stringstream gam0_ss(gam0_str);
    stringstream P0_ss(P0_str);
    stringstream seed_ss(seed_str);

    // convert to numbers
    numcells_ss >> NCELLS;
    numverts_ss >> NV;
    disp_ss >> disp;
    gam0_ss >> gam0;
    P0_ss >> P0;
    seed_ss >> seed;

    // instantiate object
    adcm2D sim(NCELLS, NV, disp, phi0, boxLengthScale, clScale, gam0, seed);

    // open output file, print
	sim.openPosObject(pos_str);

    // use contracility
    sim.useContractileShapeForce();

    // set simulation constants
    sim.setkc(kc);
    sim.setdt_pure(dt);
    sim.setgam0(gam0);
    sim.setgam(gam0);

    // compress to set pressure
    sim.nphSplitFIRE(Ftol, dPtol, P0, dt, 1);
    sim.printADCM2DConfiguration();

    // print, end
    cout << "***** NPH initial condition complete! System pressure / P0 = " << sim.getPinst() / P0 << endl;
    return 0;
}