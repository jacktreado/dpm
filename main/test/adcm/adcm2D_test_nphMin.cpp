// Test file for compression of Active Deformable Cell Model
// compile with: g++ -O3 --std=c++11 -I src main/test/adcm2D_test_nphMin.cpp  src/*.cpp -o adcm2D_test.o
// run with : ./adcm2D_test.o

#include<iostream>
#include<fstream>
#include<string>
#include<adcm2D.h>

using namespace std;

int main() {
    // input variables
    int numcells = 24;
    int numverts = 24;
    double sizeDisp = 0.1;
    double phi0 = 0.3;
    double boxLengthScale = 3.5;
    double clScale = 1.0;
    double gam0 = 1e-1;
    double kc = 1.0;
    double W = 0;
    int seed = 1;

    // NPH specific constants
    double Ftol = 1e-10;
    double P0 = 1e-4;
    double dPtol = Ftol;
    double dt0 = 0.02;

    double Posm_lambda = 0.1; // sets preferred area = a0 * (1 + Posm_lambda)
    double Posm = 1.0 / Posm_lambda;

    // output file name
    string posf = "adcm2D_test.pos";

    // initialize
    adcm2D adcm2DSimObj(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, gam0, seed);

    // open output file, print
	adcm2DSimObj.openPosObject(posf);

    // set constants
    adcm2DSimObj.setkc(kc);
    adcm2DSimObj.setW(W);
    adcm2DSimObj.setgam(gam0);
    adcm2DSimObj.setgam0(gam0);
    adcm2DSimObj.setstMat(gam0);

    // use adhesion
    double l2 = 1e-16;
    double l1 = 0.95 * l2; 
    adcm2DSimObj.useActiveTensionForce();
    adcm2DSimObj.setl1(l1);
    adcm2DSimObj.setl2(l2);

    // use osmotic pressure
    adcm2DSimObj.useOsmoticPressureShapeForce();
    adcm2DSimObj.setPosm(Posm);
    adcm2DSimObj.regularizeA0();

    // run NPH min protocol
    adcm2DSimObj.nphSplitFIRE(Ftol, dPtol, P0, dt0, true);
    // adcm2DSimObj.nphFIRE(Ftol, dPtol, P0, dt0, true);
    adcm2DSimObj.printADCM2DConfiguration();

    // To DO:
    // 1. Write new NPH-Min protocol with areal terms as CONSTRAINTS, so 
    //      -- eg: "minimize U + P_0 A subject to a_i = \bar{a}_i", where \bar_{a}_i are constraints
    //      -- Since all areas uncoupled, can just solve uncoupled quadratic equations
    // 2. Use CG instead of FIRE, I think this is cleaner
    // 3. Write simpler force and energy update functions that "shapeForces", which is a bit too convoluted at the moment
    // 4. Should be able to find configurations at a certain pressure, all with perfect areas
    // 5. Can then include box pressure into constrained area dynamics
    //      -- Why do we always need to take \lambda_+? 
    // 6. Maybe write a new code base with cell objects, easier to control where area constraints live
    //      -- Also maybe easier to implement constrained fluid pressure dynamics
    //      -- Also can push to PKS gitlab, use gitlab to handle to do list as push requests


    // run nve 
    // int NT = 2e6;
    // double T0 = 1e-6;
    // adcm2DSimObj.nve(NT, dt0, T0, true);

    return 0;
}