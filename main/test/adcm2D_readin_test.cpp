// Test file for of Active Deformable Cell Model
// compile with: g++ -O3 --std=c++11 -I src main/test/adcm2D_readin_test.cpp  src/*.cpp -o adcm2D_readin_test.o
// run with : ./adcm2D_test.o

#include<iostream>
#include<fstream>
#include<string>
#include<adcm2D.h>

using namespace std;

int main() {
    // local variables
    double dt = 0.0001;

    // input variables
    double boxLengthScale = 3.0;
    int seed = 1;
    double kc = 1.0;
    double W = 0;
    double gam = 1e-1;

    // osmotic pressure scaling
    double Posm_lambda = 0.1; // sets preferred area = a0 * (1 + Posm_lambda)
    double Posm = 1.0 / Posm_lambda;

    // input file name
    string inputf = "24cell_confluent.pos";

    // output file name
    string posf = "adcm2D_test.pos";

    // instantiate object
    adcm2D sim(inputf, seed, boxLengthScale);

    // open output file, print
	sim.openPosObject(posf);

    // set constants
    sim.setkc(kc);
    sim.setW(W);
    sim.setgam(gam);
    sim.setgam0(gam);
    sim.setstMat(gam);

    // double check pressure
    // double Ftol = 1e-8;
    // double dPtol = 1e-8;
    // double P0 = 0.5 * gam;
    // sim.nphSplitFIRE(Ftol, dPtol, P0, dt, 1);
    // return 0;

    
    // use adhesion
    double l2 = 0.1;
    double l1 = 0.5 * l2; 
    sim.useAttractiveForce();
    sim.setl1(l1);
    sim.setl2(l2);

    // use osmotic pressure
    sim.useOsmoticPressureShapeForce();
    sim.setPosm(Posm);

    // run active protocol
    double Tsim = 10;
    double Tprint = 2.0 * dt;
    double kneigh = 20.0;
    double deltaST = 0.01;

    // crawling variables
    double v0 = 0.05;
    double Ds = 0.1;
    double Dr = 0.01;

    // run
    // sim.activeBrownianCrawling(Tsim, Tprint, dt, v0, Dr, Ds);
    sim.activeTensionFluctuations(Tsim, Tprint, dt, kneigh, deltaST);

    // end
    return 0;
}