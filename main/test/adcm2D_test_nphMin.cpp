// Test file for compression of Active Deformable Cell Model
// compile with: g++ --std=c++11 -I src main/test/adcm2D_test_nphMin.cpp  src/*.cpp -o adcm2D_test.o
// run with : ./adcm2D_test.o

#include<iostream>
#include<fstream>
#include<string>
#include<adcm2D.h>

using namespace std;

int main() {
    // input variables
    int numcells = 16;
    int numverts = 24;
    double sizeDisp = 0.1;
    double phi0 = 0.5;
    double boxLengthScale = 2.75;
    double clScale = 1.0;
    double gam0 = 1e-5;
    double kc = 1.0;
    double W = 0;
    int seed = 1;

    // NPH specific constants
    double Ftol = 1e-10;
    double P0 = 2e-5;
    double dPtol = Ftol;
    double dt0 = 0.01;

    // use adhesion
    double l2 = 0.1;
    double l1 = 0.5 * l2;

    // output file name
    string posf = "adcm2D_test.pos";

    // initialize
    adcm2D adcm2DSimObj(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, gam0, seed);
    adcm2DSimObj.useAttractiveForce();

    // open output file, print
	adcm2DSimObj.openPosObject(posf);

    // set constants
    adcm2DSimObj.setkc(kc);
    adcm2DSimObj.setW(W);
    adcm2DSimObj.setgam0(gam0);
    adcm2DSimObj.setgam(gam0);

    // run NPH min protocol
    adcm2DSimObj.nphSplitFIRE(Ftol, dPtol, P0, dt0, true);
    adcm2DSimObj.printADCM2DConfiguration();

    // run nve 
    // int NT = 2e6;
    // double T0 = 1e-6;
    // adcm2DSimObj.nve(NT, dt0, T0, true);

    return 0;
}