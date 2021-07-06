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

  // specific output objects
  std::ofstream fileout;

 public:
  // constructor and destructor
  epi2D(int n, double att1, double att2, int seed)
      : dpm(n, seed) {
    z.resize(n);
    //att = attraction;
    l1 = att1;
    l2 = att2;
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
  //void setbetaEff(double val) { betaEff = val; };

  // editors & updates
  double meanl0();
  double meancalA0();
  double meankb();

  // epi cell interactions
  void vertexAttractiveForces2D();
  void attractiveForceUpdate();

  // protocols (epi needs its own protocols because of a specific use case where the sticky attraction is super strong
  // so I had to redefine attractiveForceUpdate above. only difference is the type definition of forceCall)
  //void epiCompress2Target2D(epi2DMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
  //void epiJamming2D(epi2DMemFn forceCall, double Ftol, double Ptol, double dt0, double dphi0, bool plotCompression);
  void dampedNVE2D(ofstream& enout, dpmMemFn forceCall, double B, double T, double dt0, int NT, int NPRINTSKIP);

  int getIndexOfCellLocatedHere(double xLoc, double yLoc);
  void deleteCell(double sizeRatio, int nsmall);
};

#endif