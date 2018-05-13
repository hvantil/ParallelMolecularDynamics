#include "MolDynSim.h"

MolDynSim::~MolDynSim()
{
	outProps.close();
	outPos.close();
}

double MolDynSim::PeriodicBcDist(int i,int j)
{
	double xi,yi,zi;
	atoms[i].GetPos(xi,yi,zi);
	double xj,yj,zj;
	atoms[j].GetPos(xj,yj,zj);

	Vec3 d;
	d.Set(xi-xj,yi-yj,zi-zj);
	

	if((d.x < -sz.x/2.0) && (d.x > -sz.x))
	{
		d.x += sz.x;
	}
	else if((d.x > sz.x/2.0) && (d.x < sz.x))
	{
		d.x -= sz.x;
	}

	if((d.y < -sz.y/2.0) && (d.y > -sz.y))
	{
		d.y += sz.y;
	}
	else if((d.y > sz.y/2.0) && (d.y < sz.y))
	{
		d.y -= sz.y;
	}

	if((d.z < -sz.z/2.0) && (d.z > -sz.z))
	{
		d.z += sz.z;
	}
	else if((d.z > sz.z/2.0) && (d.z < sz.z))
	{
		d.z -= sz.z;
	}

	return d.x*d.x + d.y*d.y + d.z*d.z;
}

double MolDynSim::PeriodicBcDist(int i,int j,Vec3 &d)
{
	double xi,yi,zi;
	atoms[i].GetPos(xi,yi,zi);
	double xj,yj,zj;
	atoms[j].GetPos(xj,yj,zj);

	d.Set(xi-xj,yi-yj,zi-zj);
	

	if((d.x < -sz.x/2.0) && (d.x > -sz.x))
	{
		d.x += sz.x;
	}
	else if((d.x > sz.x/2.0) && (d.x < sz.x))
	{
		d.x -= sz.x;
	}

	if((d.y < -sz.y/2.0) && (d.y > -sz.y))
	{
		d.y += sz.y;
	}
	else if((d.y > sz.y/2.0) && (d.y < sz.y))
	{
		d.y -= sz.y;
	}

	if((d.z < -sz.z/2.0) && (d.z > -sz.z))
	{
		d.z += sz.z;
	}
	else if((d.z > sz.z/2.0) && (d.z < sz.z))
	{
		d.z -= sz.z;
	}

	return d.x*d.x + d.y*d.y + d.z*d.z;
}

void MolDynSim::PeriodicBcVerlet()
{
	for(auto &a : atoms)
	{
		double rx,ry,rz;
		a.GetPos(rx,ry,rz);

		if (rx<0.0)
		{
			rx+=sz.x;
		}
		else if (rx>sz.x)
		{
			rx-=sz.x;
		}

		if (ry<0.0)
		{
			ry+=sz.y;
		}
		else if (ry>sz.y)
		{
			ry-=sz.y;
		}

		if (rz<0.0)
		{
			rz+=sz.z;
		}
		else if (rz>sz.z)
		{
			rz-=sz.z;
		}

		a.SetPos(rx,ry,rz);
	}
}


void MolDynSim::ComputeForcesSubroutine(std::vector<Vec3> &f,int jStart,int jEnd)
{
	for(int i=0;i<f.size();++i)
	{
		for(int j=jStart;j<jEnd;++j)
		{
			if(true==(i>j && (i+j)%2==0) ||
                     (i<j && (i+j)%2==1) ||
                     (i==j))
			{
				continue;
			}

			Vec3 r;
			double r2=PeriodicBcDist(i,j,r);
			if(r2<rc2)
			{
				double r1=sqrt(r2);
				double r6=r2*r2*r2;
				double r7=r6*r1;
				double r13=r6*r6*r1;
				double du_over_r = (48.0/r13 - 24.0/r7 + du_rc) / r1;
				
				double fx=du_over_r*r.x;
				double fy=du_over_r*r.y;
				double fz=du_over_r*r.z;

				f[i].x += fx;
				f[i].y += fy;
				f[i].z += fz;
				f[j].x -= fx;
				f[j].y -= fy;
				f[j].z -= fz;
			}		
		}
	}
}

void MolDynSim::ComputeForces()
{
	std::vector<int> js;
	js.push_back(0);
	for(int i=1;i<nt;++i)
	{
		js.push_back(i*(na/nt));
	}
	js.push_back(na);

	std::vector<std::vector<Vec3>> fs;
	std::vector<std::thread> threads;
	for(int i=0;i<nt;++i)
	{
		std::vector<Vec3> f;
		for(int k=0;k<na;++k)
		{
			Vec3 v;
			v.Set(0.0,0.0,0.0);
			f.push_back(v);
		}
		fs.push_back(f);
		threads.emplace_back(std::thread(
			&MolDynSim::ComputeForcesSubroutine,this,std::ref(fs[i]),js[i],js[i+1]));
	}

	for(auto &t : threads)
    {
    	t.join();
    }

    for(int k=0;k<na;++k)
    {
    	double fx=0.0;
    	double fy=0.0;
    	double fz=0.0;
    	for(int i=0;i<nt;++i)
    	{
    		fx+=fs[i][k].x;
    		fy+=fs[i][k].y;
    		fz+=fs[i][k].z;
    	}
    	atoms[k].SetForce(fx,fy,fz);
    }
}

void MolDynSim::VelocityVerlet()
{
	// update velocities (t) --> (t+dt/2)
	for(auto &a : atoms)
	{
		double vx,vy,vz,fx,fy,fz;
		a.GetVel(vx,vy,vz);
		a.GetForce(fx,fy,fz);
		a.SetVel(vx + (dt/2.0)*(fx-eta*vx),
                 vy + (dt/2.0)*(fy-eta*vy),
                 vz + (dt/2.0)*(fz-eta*vz));
	}

	// update positions (t) --> (t+dt)
	for(auto &a : atoms)
	{
		double rx,ry,rz,vx,vy,vz;
		a.GetPos(rx,ry,rz);
		a.GetVel(vx,vy,vz);
		a.SetPos(rx+dt*vx,ry+dt*vy,rz+dt*vz);
	}
	PeriodicBcVerlet();

	// update eta (t) --> (t+dt)
	ComputeKinetic();
	ComputeTemperature();
	eta += (dt/tau2) * (temp/tSet-1.0);

	// update forces (t) --> (t+dt)
	ComputeForces();

	// update velocites (t+dt/2) --> (t+dt)
	for(auto &a : atoms)
	{
		double vx,vy,vz,fx,fy,fz;
		a.GetVel(vx,vy,vz);
		a.GetForce(fx,fy,fz);
		a.SetVel((vx + dt*fx/2.0) / (1.0 + dt*eta/2.0),
                 (vy + dt*fy/2.0) / (1.0 + dt*eta/2.0),
                 (vz + dt*fz/2.0) / (1.0 + dt*eta/2.0));
	}
}


void MolDynSim::ComputeMomentum()
{
	px=0.0;
	py=0.0;
	pz=0.0;
	for(auto &a : atoms)
	{
		double x,y,z;
		a.GetVel(x,y,z);
		px+=x;
		py+=y;
		pz+=z;
	}
}

void MolDynSim::ComputePotential()
{
	pot=0.0;

	for(int i=0;i<na;++i)
	{
		for(int j=i+1;j<na;++j)
		{
			double r2=PeriodicBcDist(i,j);
			if(r2<rc2)
			{
				double r1=sqrt(r2);
				double r6=r2*r2*r2;
				double r12=r6*r6;
				double u_ij=(4.0/r12)-(4.0/r6);
				pot += u_ij - u_rc - ((r1-rc)*du_rc);
			}
		}
	}
}

void MolDynSim::ComputeKinetic()
{
	kin=0.0;
	for(auto &a : atoms)
	{
		kin+=a.GetKin();
	}
}

void MolDynSim::ComputeTemperature()
{
	// always update Kinetic before computing Temperature!
	temp=(2.0*kin)/(3.0*((double)na-1.0));
}

void MolDynSim::ComputePressure()
{
	// always update Temperature before computing Pressure!
	double fDotR=0.0;

	for(int i=0;i<na;++i)
	{
		for(int j=i+1;j<na;++j)
		{
			double r2=PeriodicBcDist(i,j);
			if(r2<rc2)
			{
				double r1=sqrt(r2);
				double r6=r2*r2*r2;
				double r7=r6*r1;
				double r13=r6*r6*r1;
				fDotR += (48.0/r13 - 24.0/r7 + du_rc) * r1;
			}
		}
	}

	pres = ((double)na*temp + fDotR/3.0) / vol;
}

void MolDynSim::WriteProps()
{
	ComputeMomentum();
	ComputePotential();
	ComputeKinetic();
	ComputeTemperature();
	ComputePressure();

	outProps << t << ",";
	outProps << px << ",";
	outProps << py << ",";
	outProps << pz << ",";
	outProps << pot << ",";
	outProps << kin << ",";
	outProps << temp << ",";
	outProps << pres << "\n";
}

void MolDynSim::WritePos()
{
	outPos << na << "\n";
	outPos << "LJ nanoparticle" << "\n";

	for(auto &a : atoms)
	{
		double x,y,z;
		a.GetPos(x,y,z);
		outPos << "H " << x << " " << y << " " << z << "\n";
	}
}


void MolDynSim::SetParams(int nIters,int nThreads,int outItvl,double dTime)
{
	ni=nIters;
	nt=nThreads;
	oi=outItvl;
	dt=dTime;
}

void MolDynSim::SetProps(int nAtoms,double rCutoff,double tempSet,double timeConst,
                         double boxSizeX,double boxSizeY,double boxSizeZ)
{
	na=nAtoms;
	rc=rCutoff;
	tSet=tempSet;
	tau=timeConst;

	sz.Set(boxSizeX,boxSizeY,boxSizeZ);

	for(int i=0;i<na;++i)
	{
		Atom a;
		atoms.push_back(a);
	}

	rc2=rc*rc;
	tau2=tau*tau;
	du_rc=(24.0/pow(rc,7))-(48.0/pow(rc,13));
	u_rc=(4.0/pow(rc,12))-(4.0/pow(rc,6));
	vol=sz.x*sz.y*sz.z;

	t=0.0;
	eta=0.0;
	px=0.0;
	py=0.0;
	pz=0.0;
	pot=0.0;
	kin=0.0;
	temp=0.0;
	pres=0.0;
}

void MolDynSim::Initialize()
{
	// setup output files
	std::string fn1="outProps";
	fn1.append(std::to_string(na));
	fn1.append(".csv");
	outProps.open(fn1.data());

	std::string fn2="outPos";
	fn2.append(std::to_string(na));
	fn2.append(".xyz");
	outPos.open(fn2.data());

	// read positions from input file
	std::string fn3="in";
	fn3.append(std::to_string(na));
	fn3.append(".txt");

	std::ifstream in;
	in.open(fn3.data());
	for(auto &a : atoms)
	{
		double x,y,z;
		in >> x;
		in >> y;
		in >> z;
		a.SetPos(x,y,z);
	}
	in.close();

	// random inital velocities [-1, 1]
	for(auto &a : atoms)
	{
		a.SetVel(2.0*(((double)rand()/(double)RAND_MAX)-0.5),
                 2.0*(((double)rand()/(double)RAND_MAX)-0.5),
                 2.0*(((double)rand()/(double)RAND_MAX)-0.5));
	}

	// shift to get zero total momentum
	ComputeMomentum();
	double pxAvg=px/(double)na;
	double pyAvg=py/(double)na;
	double pzAvg=pz/(double)na;
	for(auto &a : atoms)
	{
		double x,y,z;
		a.GetVel(x,y,z);
		a.SetVel(x-pxAvg,y-pyAvg,z-pzAvg);
	}

	// scale to get temp = tempSet	
	ComputeKinetic();
	ComputeTemperature();	
	double alpha=sqrt(tSet/temp);	
	for(auto &a : atoms)
	{
		double x,y,z;
		a.GetVel(x,y,z);
		a.SetVel(x*alpha,y*alpha,z*alpha);
	}

	// set forces to zero
	for(auto &a : atoms)
	{
		a.SetForce(0.0,0.0,0.0);
	}
}

void MolDynSim::Run()
{
	for(int i=0;i<ni;++i)
	{
		if(i%(ni/25) == 0)
		{
			printf("time step: %d\n",i);
		}

		if(i%oi == 0)
		{
			t=i*dt;
			WritePos();
			WriteProps();
		}

		VelocityVerlet();
	}
}

