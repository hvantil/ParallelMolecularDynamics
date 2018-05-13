#ifndef ATOM_H
#define ATOM_H

#include "Vec3.h";
#include <stdio.h>

class Atom
{
private:
	Vec3 r,v,f;

public:
	Atom();
	void SetPos(double x,double y,double z);
	void SetVel(double x,double y,double z);
	void SetForce(double x,double y,double z);
	void GetPos(double &x,double &y,double &z);
	void GetVel(double &x,double &y,double &z);
	void GetForce(double &x,double &y,double &z);
	double GetKin();
};

#endif