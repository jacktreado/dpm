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
    int numcells = 12;
    int numverts = 16;
    double sizeDisp = 0.1;
    double phi0 = 0.45;
    double boxLengthScale = 3.0;
    double clScale = 0.5;
    double gam0 = 0.01;
    double kc = 0.5;
    int seed = 1;

    // NPH specific constants
    double Ftol = 1e-12;
    double P0 = 1e-6;
    double dPtol = Ftol;
    double dt0 = 0.01;

    // output file name
    string posf = "adcm2D_test.pos";

    // initialize
    adcm2D adcm2DSimObj(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, gam0, seed);

    // open output file, print
	adcm2DSimObj.openPosObject(posf);

    // set constants
    adcm2DSimObj.setkc(kc);

    // run NPH min protocol
    adcm2DSimObj.nphFIRE(Ftol, dPtol, P0, dt0, true);
    adcm2DSimObj.printADCM2DConfiguration();

    // run nve 
    // int NT = 2e6;
    // double T0 = 1e-6;
    // adcm2DSimObj.nve(NT, dt0, T0, true);

    return 0;
}