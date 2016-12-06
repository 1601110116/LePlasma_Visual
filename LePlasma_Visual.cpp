//============================================================================
// Name        : ZLOMP.cpp
// Author      : ZL
// Version     :
// Copyright   : ZL-CSPIC
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <CaseLauncher.h>
#include <DispersionRelation.h>
#include <UniformE.h>
#include <UniformB.h>
#include <UniformBFalse.h>
#include <GLEnvironment.h>
#include <Grid.h>
#include <LePlasma.h>
#include <RunManager.h>
#include <Visualize.h>


using namespace std;

int main(int argc, char**argv){


#if MPI_PARALLEL
	MPI_Init(&argc, &argv);

	omp_set_num_threads(1);

	RunManager::MPI_Init_Prarllel();
#endif


	Case* myCase=new DispersionRelation();
	//Case* myCase=new SingleParticle();

	RunManager *run = new RunManager(myCase);

	InitGLEnvironment();

#if MPI_PARALLEL
	Visualize* stage=new Visualize(myCase->grid,RunManager::MPI_ID);
#else
	Visualize* stage=new Visualize(myCase->grid);
#endif

	stage->idleFn=run->stepNext;

	InitStage(stage);


#if MPI_PARALLEL
	MPI_Finalize();
#endif

	return 0;
}
