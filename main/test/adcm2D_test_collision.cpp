// Test file for collision of 2 Active Deformable Cells (passive, just repulsive vertices)
// compile with: g++ --std=c++11 -I src main/test/adcm2D_test_collision.cpp  src/*.cpp -o adcm2D_test_collision.o
// run with : ./adcm2D_test_collision.o

// TO DO:
// * Check energy conservation with normal area term (because there is a well defined energy there) OR add appropriate energy for osmotic pressure force

#include<iostream>
#include<fstream>
#include<string>
#include<adcm2D.h>

using namespace std;

int main() {
    // input variables
    int numcells = 3;
    int numverts = 32;
    double sizeDisp = 0;
    double phi0 = 0.2;
    double boxLengthScale = 3.5;
    double clScale = 1.0;
    double gam0 = 0.1;
    double kc = 0.1;
    int seed = 2;

    // collision variables
    double Tend = 500;
    double Tskip = 1.0;
    double dt = 0.001;
    double v0 = 0.05;
    double cx, cy, vx, vy;    

    // output file name
    string posf = "adcm2D_test_collision.pos";
    string enf = "adcm2D_test_collision.en";
    string forcef = "adcm2D_test_collision.frc";

    // open force output file
    ofstream forceOuputObj(forcef.c_str());

    // initialize
    adcm2D sim(numcells, numverts, sizeDisp, phi0, boxLengthScale, clScale, gam0, seed);

    // open output file, print
	sim.openPosObject(posf);

    // set constants
    sim.setkc(kc);
    sim.setdt_pure(dt);
    sim.setgam(gam0);
    sim.setgam0(gam0);
    sim.setstMat(gam0);

     // use adhesion
    double l2 = 0.05;
    double l1 = 0.5 * l2; 
    sim.useActiveTensionForce();
    sim.setl1(l1);
    sim.setl2(l2);
    // sim.useRepulsiveForce();

    // use osmotic pressure
    sim.useOsmoticPressureShapeForce();
    sim.setPosm(10.0);
    sim.regularizeA0();

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
    double t = 0.0;
    int nprint = 0;
    double Px, Py;
    while (t < Tend){
        // VV velocity half-step update
		sim.vvVelUpdate();

		// VV positions full-step update
		sim.vvPosUpdate();

		// Update Forces & Potential Energy
		sim.activeTensionForceUpdate();
        // sim.adcm2DForceUpdate();

		// VV velocity half-step update
		sim.vvVelUpdate();

        // print statement
		if (floor(t / Tskip) == nprint){
            // get current linear momentum
            sim.netLinearMomentum2D(Px, Py);

			cout << endl << endl;
			cout << "===========================================" << endl;
			cout << " 	N V E 									" << endl;
			cout << "===========================================" << endl;
			cout << endl;
            cout << "   ** t        = " << t << endl;
			cout << "   ** nprint 	= " << nprint << endl;
			cout << "	** U 		= " << sim.getU() << endl;
			cout << " 	** K 		= " << sim.vertexKineticEnergy() << endl;
			cout << " 	** E 		= " << sim.getU() + sim.vertexKineticEnergy() << endl;
            cout << "   ** net Px   = " << Px << endl;
            cout << "   ** net Py   = " << Py << endl;
            cout << "   ** net Lz   = " << sim.netAngularMomentum2D() << endl << endl;

			// print it
            sim.printADCM2DConfiguration();
            sim.printInstantaneousForces(forceOuputObj);
            nprint++;
		}

        // update time
        t += dt;
    }

    forceOuputObj.close();
    return 0;
}