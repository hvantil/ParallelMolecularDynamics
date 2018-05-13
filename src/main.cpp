#include <stdlib.h>
#include <time.h>
#include "fssimplewindow.h"
#include "MolDynSim.h"

int main(void)
{
	srand(time(NULL));
	FsChangeToProgramDir();

	MolDynSim sim;

	int nIters=10000;
	int nThreads=1;
	int outItvl=10;
	double dTime=0.002;
	sim.SetParams(nIters,nThreads,outItvl,dTime);

	int nAtoms=256;
	double rCutoff=2.5;
	double tempSet=100.0/120.96;
	double timeConst=0.05;
	double boxSizeX=6.8;
	double boxSizeY=6.8;
	double boxSizeZ=6.8;
	sim.SetProps(nAtoms,rCutoff,tempSet,timeConst,
                 boxSizeX,boxSizeY,boxSizeZ);

	sim.Initialize();
	
	clock_t tStart=clock();
	sim.Run();
	clock_t tEnd=clock();
	printf("runtime = %f sec\n", (tEnd-tStart)/(double)CLOCKS_PER_SEC);

	return 0;
}
