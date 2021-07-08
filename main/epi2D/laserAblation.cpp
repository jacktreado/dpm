// File to compress DPM into confluent epithelial layer, then laser ablate
// ** Features to add: crawling, purse-string contraction, substrate interaction
//
// Cells are bidisperse, with bending energy (broken) and vertex-vertex attraction options.
//
//
// Create bidisperse DPM particles, set constants, place particle centers,
// relax shapes + positions, compress to target phi, delete cells, conduct
// damped MD dynamics.
//
// Compilation command:
// g++ -O3 --std=c++11 -I src main/epi2D/laserAblation.cpp src/dpm.cpp src/epi2D.cpp -o main/epi2D/laserAblation.o
// ./main/epi2D/laserAblation.o 12 20 1.08 0.7 0.9 1.0 0.5 0.5 1.0 1 100 pos.test energy.test
//
//
// Parameter input list
// 1. NCELLS: 			number of particles
// 2. nsmall: 			number of vertices on small particles (larger
// particles set by 1.4:1.0 bidispersity)
// 3. calA0: 			  preferred shape parameter for all particles
// 4. phiMin: 			p ?
// 5. phiMax: 			p
// 6. kl: 				  perimeter spring constant
// 7. att:          attraction range and strength parameter
// 8. B:            damping coefficient
// 9. Dr0:          rotational diffusion constant for protrusion activity
// 10. seed: 			  seed for random number generator
// 11. time:        amount of time (tau) to simulate
// 12. positionFile: 	string of path to output file with position/configuration data
// 13. energyFile:  string of path to output file with energy data

// header files
#include <sstream>
#include "epi2D.h"

// preprocessor macros
#define NDIM 2

// namspace
using namespace std;

const bool plotCompression = 0;     // whether or not to plot configuration during compression protocol
const double dphi0 = 0.005;         // packing fraction increment
const double ka = 1.0;              // area force spring constant (should be unit)
const double kc = 1.0;              // interaction force spring constant (should be unit)
const double boxLengthScale = 2.5;  // neighbor list box size in units of initial l0
//const double phi0 = 0.5;            // initial packing fraction
const double smallfrac = 0.5;  // fraction of small particles
const double sizeratio = 1.4;  // size ratio between small and large particles
const double dt0 = 0.005;      // initial magnitude of time step in units of MD time
const double Ptol = 1e-8;
const double Ftol = 1e-12;

int main(int argc, char const* argv[]) {
  // local variables to be read in
  int NCELLS, nsmall, seed;
  double calA0, kl, kb, phiMin, phiMax, att, B, NT_dbl, Dr0, time_dbl;

  // read in parameters from command line input
  string NCELLS_str = argv[1];
  string nsmall_str = argv[2];
  string calA0_str = argv[3];
  string phiMin_str = argv[4];
  string phiMax_str = argv[5];
  string kl_str = argv[6];
  string att_str = argv[7];
  string B_str = argv[8];
  string Dr0_str = argv[9];
  string seed_str = argv[10];
  string time_str = argv[11];
  string positionFile = argv[12];
  string energyFile = argv[13];

  // using sstreams to get parameters
  stringstream NCELLSss(NCELLS_str);
  stringstream nsmallss(nsmall_str);
  stringstream calA0ss(calA0_str);
  stringstream phiMinss(phiMin_str);
  stringstream phiMaxss(phiMax_str);
  stringstream klss(kl_str);
  stringstream Dr0ss(Dr0_str);
  stringstream attss(att_str);
  stringstream Bss(B_str);
  stringstream seedss(seed_str);
  stringstream timess(time_str);

  // read into data
  NCELLSss >> NCELLS;
  nsmallss >> nsmall;
  calA0ss >> calA0;
  phiMinss >> phiMin;
  phiMaxss >> phiMax;
  klss >> kl;
  Dr0ss >> Dr0;
  attss >> att;
  Bss >> B;
  seedss >> seed;
  timess >> time_dbl;

  std::ofstream enout(energyFile);

  // number of time steps
  int NT = int(time_dbl / dt0);

  const double phi0 = phiMin;

  // instantiate object
  epi2D epithelial(NCELLS, 0.0, 0.0, Dr0, seed);
  epithelial.setl1(att);
  epithelial.setl2(att / 2);

  epithelial.openPosObject(positionFile);

  // set spring constants
  epithelial.setka(ka);
  epithelial.setkl(kl);
  epithelial.setkb(kb);
  epithelial.setkc(kc);

  epithelial.bidisperse2D(calA0, nsmall, smallfrac, sizeratio);

  cout << "initializePositions2D\n";
  epithelial.initializePositions2D(phi0, Ftol);

  cout << "initializeNeighborLinkedList2D\n";
  epithelial.initializeNeighborLinkedList2D(boxLengthScale);

  // set base dpm force, upcast derived epi2D forces
  dpmMemFn repulsiveForceUpdate = &dpm::repulsiveForceUpdate;
  dpmMemFn attractiveForceUpdate = static_cast<void (dpm::*)()>(&epi2D::attractiveForceUpdate_2);
  dpmMemFn activeForceUpdate = static_cast<void (dpm::*)()>(&epi2D::activeAttractiveForceUpdate);

  cout << "starting vertexCompress2Target2D\n";
  epithelial.vertexCompress2Target2D(repulsiveForceUpdate, Ftol, dt0, phiMax, dphi0);

  //after compress, turn on damped NVE
  double T = 1e-4;
  epithelial.drawVelocities2D(T);
  epithelial.dampedNVE2D(enout, attractiveForceUpdate, B, dt0, NT / 10, NT / 10);

  /*********** placeholder: LASER ABLATION SCHEME ************************/
  double xLoc = 0.0, yLoc = 0.0;
  int numCellsToAblate = 3;
  //epithelial.laserAblate(numCellsToAblate, sizeratio, nsmall, xLoc, yLoc);
  cout << "numCells = " << epithelial.getNCELLS() << '\n';
  for (int ci = 0; ci < epithelial.getNCELLS(); ci++) {
    cout << "Psi(" << ci << ") = " << epithelial.getPsi(ci) << '\n';
    //epithelial.orientDirector(ci, 0.0, 0.0);
    cout << "Psi(" << ci << ") = " << epithelial.getPsi(ci) << '\n';
  }
  epithelial.dampedNVE2D(enout, activeForceUpdate, B, dt0, NT, NT / 10);
  for (int ci = 0; ci < epithelial.getNCELLS(); ci++) {
    cout << "Psi(" << ci << ") = " << epithelial.getPsi(ci) << '\n';
    //epithelial.orientDirector(ci, 0.0, 0.0);
    cout << "Psi(" << ci << ") = " << epithelial.getPsi(ci) << '\n';
  }

  /*********placeholder: ELASTIC SHEET FRACTURE SCHEME*********************/
  //epithelial.holePunching(sizeratio, nsmall, enout, attractiveForceUpdate, B, =1000, 1000);

  cout << "\n** Finished laserAblation.cpp, ending. " << endl;

  //testing energy conservation
  enout.close();

  return 0;
}