#ifndef MOL_DYN_SIM_H
#define MOL_DYN_SIM_H

#include <vector>
#include <string>
#include <fstream>
#include <thread>
#include <functional> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Vec3.h"
#include "Atom.h"

class MolDynSim
{
private:
	int ni,nt,oi;
	double dt;

	int na;
	double rc,tSet,tau;
	double rc2,tau2,du_rc,u_rc,vol;
	double t,eta,px,py,pz,pot,kin,temp,pres;

	Vec3 sz;
	std::vector<Atom> atoms;

	std::ofstream outProps, outPos;

private:
	double PeriodicBcDist(int i,int j);
	double PeriodicBcDist(int i,int j,Vec3 &d);
	void PeriodicBcVerlet();

	void ComputeForcesSubroutine(std::vector<Vec3> &f,int jStart,int jEnd);
	void ComputeForces();
	void VelocityVerlet();

	void ComputeMomentum();
	void ComputePotential();
	void ComputeKinetic();
	void ComputeTemperature();
	void ComputePressure();
	void WriteProps();
	void WritePos();

public:
	~MolDynSim();
	void SetParams(int nIters,int nThreads,int outItvl,double dTime);
	void SetProps(int nAtoms,double rCutoff,double tempSet,double timeConst,
                  double boxSizeX,double boxSizeY,double boxSizeZ);
	void Initialize();
	void Run();
};

#endif