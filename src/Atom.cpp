#include "Atom.h"

Atom::Atom()
{
	r.Set(0.0,0.0,0.0);
	v.Set(0.0,0.0,0.0);
	f.Set(0.0,0.0,0.0);
}

void Atom::SetPos(double x,double y,double z)
{
	r.Set(x,y,z);
}

void Atom::SetVel(double x,double y,double z)
{
	v.Set(x,y,z);
}

void Atom::SetForce(double x,double y,double z)
{
	f.Set(x,y,z);
}

void Atom::GetPos(double &x,double &y,double &z)
{
	x=r.x;
	y=r.y;
	z=r.z;
}

void Atom::GetVel(double &x,double &y,double &z)
{
	x=v.x;
	y=v.y;
	z=v.z;
}

void Atom::GetForce(double &x,double &y,double &z)
{
	x=f.x;
	y=f.y;
	z=f.z;
}

double Atom::GetKin()
{
	return (v.x*v.x + v.y*v.y + v.z*v.z)/2.0;
}
