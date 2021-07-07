#ifndef EPI2D_H
#define EPI2D_H

/*

	HEADER FILE FOR epi CLASS

		-- Inherits from DPM class 
		-- For collections of epithelial cells
		-- Incorporates ?
		-- ONLY FOR 2D

*/

#include <algorithm>
#include <stdexcept>
#include "dpm.h"

// namespace
using namespace std;

class epi2D;
typedef void (epi2D::*epi2DMemFn)(void);

class epi2D : public dpm {
 protected:
  // bending energy per vertex
  // NOTE: will need to add different Hessian computation
  std::vector<double> kbi;

  // vertex-vertex contact network
  std::vector<bool> gij;

  // vv contacts per cell
  std::vector<int> z;

  // adhesion strength
  //double att;

  // directors (direction vectors) for activity
  std::vector<double> psi;

  // rotational diffusion constant
  double Dr0;

  // specific output objects
  std::ofstream fileout;



 public:
  // constructor and destructor
  epi2D(int n, double att1, double att2, double Dr, int seed)
      : dpm(n, seed) {
    z.resize(n);
    //att = attraction;
    l1 = att1;
    l2 = att2;
    Dr0 = Dr;
    vector<double> temp(NCELLS, 0.0);
    psi = temp;
    for (int ci = 0; ci < NCELLS; ci++)
      psi.at(ci) = 2 * PI * (drand48() - 0.5);
  };

  // File openers
  void openFileObject(std::string& str) {
    fileout.open(str.c_str());
    if (!fileout.is_open()) {
      std::cout << "	ERROR: file could not open " << str << "..." << std::endl;
      exit(1);
    } else
      std::cout << "** Opening file " << str << " ..." << std::endl;
  }

  // setters
  void setkbi(double val) { fill(kbi.begin(), kbi.end(), val); };

  // editors & updates
  double meanl0();
  double meancalA0();
  double meankb();

  // epi cell interactions
  void vertexAttractiveForces2D_2();
  void attractiveForceUpdate_2();
  void activeAttractiveForceUpdate();

  // protocols
  void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  void zeroMomentum();
  void dampedNVE2D(ofstream& enout, dpmMemFn forceCall, double B, double dt0, int NT, int NPRINTSKIP);

  int getIndexOfCellLocatedHere(double xLoc, double yLoc);
  void deleteCell(double sizeRatio, int nsmall, double xLoc, double yLoc);
  void laserAblate(int numCellsAblated, double sizeRatio, int nsmall, double xLoc, double yLoc);

  void isotropicDistanceScaling(ofstream& enout, dpmMemFn forceCall, double B, int NT, int NPRINTSKIP);
  void holePunching(double sizeRatio, int nsmall, ofstream& enout, dpmMemFn forceCall, double B, int NT, int NPRINTSKIP);
};

#endif