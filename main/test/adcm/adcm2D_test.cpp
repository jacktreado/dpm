// Test file for of Active Deformable Cell Model
// compile with: g++ --std=c++11 -I src main/test/adcm2D_test.cpp  src/*.cpp -o adcm2D_test.o
// run with : ./adcm2D_test.o

#include<iostream>
#include<fstream>
#include<string>
#include<adcm2D.h>

using namespace std;

int main() {
    // input variables
    int numcells = 64;
    int numverts = 16;
    double sizeDisp = 0.1;
    double phi0 = 0.8;
    double boxLengthScale = 3.0;
    double clScale = 1.0;
    double gam0 = 0.01;
    int seed = 1;
    double kc = 1.0;
    double W = 0.0;

    // output file name
    string posf = "adcm2D_test.pos";

    // instantiate object
    adcm2D sim(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, gam0, seed);

    // open output file, print
	sim.openPosObject(posf);

    // set constants
    sim.setkc(kc);
    sim.setW(W);

    // use adhesion
    double l2 = 0.02;
    double l1 = 0.5 * l2;
    sim.useAttractiveForce();
    sim.setl1(l1);
    sim.setl2(l2);

    // run ABC protocol
    double Tsim = 2000.0;
    double Tprint = 5.0;
    double v0 = 0.01;
    double Ds = 0.1;
    double Dr = 0.01;
    double dt = 0.01;

    // run
    sim.activeBrownianCrawling(Tsim, Tprint, dt, v0, Dr, Ds);

    // end
    return 0;
}