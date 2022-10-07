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
    int numcells = 10;
    int numverts = 24;
    double sizeDisp = 0.1;
    double phi0 = 0.3;
    double boxLengthScale = 3.0;
    double clScale = 0.1;
    double gam0 = 0.1;
    int seed = 1;

    // output file name
    string posf = "adcm2D_test.pos";

    // initialize
    adcm2D adcm2DSimObj(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, gam0, seed);

    // open output file, print
	adcm2DSimObj.openPosObject(posf);
    adcm2DSimObj.printADCM2DConfiguration();

    return 0;
}