#ifndef DPM_H
#define DPM_H

/*

	HEADER FILE FOR DPM CLASS

		-- Stores vertex information for NCELLS with n_\mu vertices each
		-- Stores shape parameters + spring constants
		-- Integrates equations of motion (NVE/NVT)
		-- Finds specified configurations (jammed, fixed p, fixed shear strain, etc)
		-- Compute mechanical information (Hessian, elastic moduli, etc)
		-- Prints data to files

	Jack Treado, 04/10/21

*/


#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <functional>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// pointer-to-member function call macro
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

class dpm;
typedef void (dpm::*dpmMemFn)(void);

// global constants
const double PI = 4.0 * atan(1.0);
const int nvmin = 12;

// printing constants
const int w = 10;
const int wnum = 25;
const int pnum = 14;

// FIRE constants
const double alpha0 = 0.2;
const double finc = 1.1;
const double fdec = 0.5;
const double falpha = 0.99;

const int NSKIP = 20000;
const int NMIN = 10;
const int NNEGMAX = 1000;
const int NDELAY = 20;
const int itmax = 5e7;

class dpm
{
protected:
	// int scalars
	int NCELLS;
	int NDIM;
	int NNN;
	int NVTOT;
	int vertDOF;

	// time step size
	double dt;

	// potential energy
	double U;

	// particle spring constants
	double ka;
	double kl;
	double kb;
	double kc;

	// particle attraction constants
	double l1, l2;

	// boundary parameters
	std::vector<double> L;
	std::vector<bool> pbc;

	// particle shape parameters
	std::vector<double> a0;
	std::vector<double> l0;
	std::vector<double> t0;
	std::vector<double> r;

	// indexing variables
	std::vector<int> nv;
	std::vector<int> szList;
	std::vector<int> im1;
	std::vector<int> ip1;

	// dynamical variables
	std::vector<double> x;
	std::vector<double> v;
	std::vector<double> F;

	// macroscopic stress vector
	std::vector<double> stress;

	// contact network
	std::vector<int> cij;

	// Box linked-list variables
	int NBX;
	std::vector<int> sb;
	std::vector<double> lb;
	std::vector<std::vector<int>> nn;
	std::vector<int> head;
	std::vector<int> last;
	std::vector<int> list;

	// output objects
	std::ofstream posout;

public:
	// Constructors and Destructors
	dpm(int ndim);
	dpm(int n, int ndim, int seed);
	dpm(int n, int seed) : dpm(n, 2, seed) {}
	~dpm();

	// -- G E T T E R S

	// main ints
	int getNCELLS() { return NCELLS; };
	int getNDIM() { return NDIM; };
	int getNNN() { return NNN; };
	int getNVTOT() { return NVTOT; };
	int getvertDOF() { return vertDOF; };
	int getNV(int ci) { return nv.at(ci); };

	// force parameters
	double getdt() { return dt; };
	double getka() { return ka; };
	double getkl() { return kl; };
	double getkb() { return kb; };
	double getkc() { return kc; };

	// static cell info
	double geta0(int ci) { return a0[ci]; };
	double getl0(int gi) { return l0[gi]; };
	double gett0(int gi) { return t0[gi]; };
	double getr(int gi) { return r[gi]; };

	// dynamic cell info
	double getx(int gi, int d) { return x[NDIM * gi + d]; };
	double getv(int gi, int d) { return v[NDIM * gi + d]; };
	double getF(int gi, int d) { return F[NDIM * gi + d]; };
	double getU() { return U; };

	// boundary variables
	double getL(int d) { return L.at(d); };
	bool getpbc(int d) { return pbc.at(d); };

	// cell shape indexing + information
	int gindex(int ci, int vi);
	void cindices(int &ci, int &vi, int gi);
	double area(int ci);
	double perimeter(int ci);
	void com2D(int ci, double &cx, double &cy);
	double vertexPackingFraction2D();
	double vertexPreferredPackingFraction2D();
	double vertexKineticEnergy();
	int vvContacts();
	int ccContacts();

	// Setters
	void setpbc(int d, bool val) { pbc.at(d) = val; };
	void setNCELLS(int val) { NCELLS = val; };
	void setdt(double val);
	void setka(double val) { ka = val; };
	void setkl(double val) { kl = val; };
	void setkb(double val) { kb = val; };
	void setkc(double val) { kc = val; };
	void setl1(double val) { l1 = val; };
	void setl2(double val) { l2 = val; };

	// File openers
	void openPosObject(std::string &str)
	{
		posout.open(str.c_str());
		if (!posout.is_open()) {
			std::cerr << "	ERROR: posout could not open " << str << "..." << std::endl;
			exit(1);
		}
		else
			std::cout << "** Opening pos file " << str << " ..." << std::endl;
	}
    
	// Initialize particles (two dimensions)
	void monodisperse2D(int n);
	void bidisperse2D(double calA0, int nsmall, double smallfrac, double sizefrac);
	void gaussian2D(double dispersion, double calA0, int n1);
	void sinusoidalPreferredAngle(double thA, double thK);
	void initializeVertexShapeParameters(double calA0, int nref);
	void initializeVertexShapeParameters(int ci, double calA0, double lenscale);
	void initializeVertexIndexing2D();
	void initializePositions2D(double phi0, double Ftol);
	void initializeNeighborLinkedList2D(double boxLengthScale);

	// editing & updating
	void sortNeighborLinkedList2D();
	void scaleParticleSizes2D(double scaleFactor);
	int removeRattlers();
	void drawVelocities2D(double T);

	// force definitions
	void resetForcesAndEnergy();
	void shapeForces2D();
	void vertexRepulsiveForces2D();
	void vertexAttractiveForces2D();

	// force updates
	void repulsiveForceUpdate();
	void attractiveForceUpdate();

	// simple integrators
	void vertexFIRE2D(dpmMemFn forceCall, double Ftol, double dt0);
	void vertexNVE2D(std::ofstream &enout, dpmMemFn forceCall, double T, double dt0, int NT, int NPRINTSKIP);
	void vertexLangevinNVT2D(std::ofstream &enout, dpmMemFn forceCall, double T0, double gam, double dt0, int NT, int NPRINTSKIP);
	void vertexLangevinNVT2D(dpmMemFn forceCall, double T0, double gam, double dt0, int NT, int NPRINTSKIP);

	// protocols
	void vertexCompress2Target2D(dpmMemFn forceCall, double Ftol, double dt0, double phi0Target, double dphi0);
	void vertexJamming2D(dpmMemFn forceCall, double Ftol, double Ptol, double dt0, double dphi0, bool plotCompression);
	void vertexAnneal2Jam2D(dpmMemFn forceCall, double Ftol, double Ptol, double dt0, double dphi0, double T0, double trun, bool plotCompression);

	// hessian methods
	// note: dynamical matrix contribution is always M = H - S
	void dpmHessian2D(Eigen::MatrixXd &H, Eigen::MatrixXd &S);
	void dpmAreaHessian2D(Eigen::MatrixXd &Ha, Eigen::MatrixXd &Sa);
	void dpmPerimeterHessian2D(Eigen::MatrixXd &Hl, Eigen::MatrixXd &Sl);
	void dpmBendingHessian2D(Eigen::MatrixXd &Hb, Eigen::MatrixXd &Sb);
	void dpmRepulsiveHarmonicSprings2D(Eigen::MatrixXd &Hvv, Eigen::MatrixXd &Svv);

	// print vertex information to file
	void printNeighborList();
	void printContactMatrix();
	void printConfiguration2D();
	void printHessianEigenvalues2D(std::ofstream &hessout, Eigen::MatrixXd &M);
};

#endif
