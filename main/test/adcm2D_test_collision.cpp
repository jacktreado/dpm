// Test file for collision of 2 Active Deformable Cells (passive, just repulsive vertices)
// compile with: g++ --std=c++11 -I src main/test/adcm2D_test_collision.cpp  src/*.cpp -o adcm2D_test_collision.o
// run with : ./adcm2D_test_collision.o

#include<iostream>
#include<fstream>
#include<string>
#include<adcm2D.h>

using namespace std;

int main() {
    // input variables
    int numcells = 4;
    int numverts = 16;
    double sizeDisp = 0;
    double phi0 = 0.2;
    double boxLengthScale = 3.0;
    double clScale = 0.5;
    double gam0 = 0.01;
    double kc = 0.5;
    int seed = 10;

    // collision variables
    int NT = 2e5;
    int NSKIP = 2e3;
    double dt = 0.005;
    double v0 = 0.1;
    double cx, cy, vx, vy;

    // output file name
    string posf = "adcm2D_test_collision.pos";
    string enf = "adcm2D_test_collision.en";

    // initialize
    adcm2D sim(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, gam0, seed);

    // open output file, print
	sim.openPosObject(posf);

    // set constants
    sim.setkc(kc);
    sim.setdt_pure(dt);

    // initial collision velocities
    for (int ci=0; ci<numcells; ci++){
        // get center of mass
        sim.com2D(ci, cx, cy);

        // get vector to center of box
        vx = 0.5 * sim.getL(0) - cx;
        vy = 0.5 * sim.getL(1) - cy;
        
        // rescale
        vx = vx / sqrt(vx * vx + vy * vy);
        vy = vy / sqrt(vx * vx + vy * vy);
        vx = vx * v0;
        vy = vy * v0;

        // set initial velocities
        sim.setComV(ci, 0, vx);
        sim.setComV(ci, 1, vy);
    }

    // sim.printNeighborList();

    // integrate time
    for (int tt=0; tt<NT; tt++){
        // VV velocity half-step update
		sim.vvVelUpdate();

		// VV positions full-step update
		sim.vvPosUpdate();

		// Update Forces & Potential Energy
		sim.adcm2DForceUpdate();

		// VV velocity half-step update
		sim.vvVelUpdate();

        // print statement
		if (tt % NSKIP == 0){
			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << " 	N V E 									" << endl;
			cout << "===========================================" << endl;
			cout << endl;
			cout << "	** tt 		= " << tt << endl;
			cout << "	** U 		= " << sim.getU() << endl;
			cout << " 	** K 		= " << sim.vertexKineticEnergy() << endl;
			cout << " 	** E 		= " << sim.getU() + sim.vertexKineticEnergy() << endl;
			cout << "	** dt 		= " << dt << endl;

			// print it
            sim.printADCM2DConfiguration();
		}
    }

    return 0;
}