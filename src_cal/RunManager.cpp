/*
 * RunManager.cpp
 *
 *  Created on: Mar 6, 2016
 *      Author: zlstudio
 */

#include <Case.h>
#include <Engine.h>
#include <LePlasma.h>
#include <mpi.h>
#include <MPIGrid.h>
#include <omp.h>
#include <Range.h>
#include <RunManager.h>
#include <iostream>

Grid* RunManager::grid=NULL;
double RunManager::deltaT=0.0;
void (*RunManager::Method)(const Range&, double)=NULL;
int RunManager::MPI_ID=0;
Range* RunManager::ranges;
int RunManager::Nodes;
Case* RunManager::_Case;

static int n=0;
static double dt1=0;
static double dt2=0;
static double dt3=0;

using namespace std;


RunManager::RunManager(MPIGrid* _grid,void (*_Method)(const Range &,double),double _dt) {
	// TODO Auto-generated constructor stub
	grid=_grid;
	Method=_Method;
	deltaT=_dt;

	Init_Parallel();
}

RunManager::RunManager(Case* _case){
	grid=_case->grid;
	deltaT=_case->deltaT;
	_Case=_case;

	Init_Parallel();
}
RunManager::~RunManager() {
	// TODO Auto-generated destructor stub
}

void RunManager::stepNext(){

	n++;

	double t=omp_get_wtime();

#pragma omp parallel
	{
		if(_Case!=NULL){
			_Case->engine->update(ranges[omp_get_thread_num()]);
		}else{
			Method(ranges[omp_get_thread_num()],deltaT);
		}
	}

	dt1+=(omp_get_wtime()-t)*1000;

	t=omp_get_wtime();

	grid->refreshParticleLocation();

	dt2+=(omp_get_wtime()-t)*1000;

#if REPORT

	t=omp_get_wtime();

	grid->reportToHost();

	dt3+=(omp_get_wtime()-t)*1000;

#endif

#if RECORD
	if(_Case!=NULL){
		_Case->report();
	}

#endif

	if(n>=10){
		cout<<"ID: "<<MPI_ID<<" Time:\tCal\t"<<dt1/n<<" ms\tArr:\t"<<dt2/n<<" ms\t Report:\t"<<dt3/n<<endl;
		dt1=dt2=dt3=n=0;
	}
}

/*****INIT MPI*****/
void RunManager::MPI_Init_Prarllel(){

	int l=30;
	char name[30];
	MPI_Comm_rank(MPI_COMM_WORLD,&MPI_ID);
	MPI_Comm_size(MPI_COMM_WORLD,&Nodes);
	MPI_Get_processor_name(name,&l);

	int prevNode=(MPI_ID==0?Nodes-1:MPI_ID-1);
	int nextNode=(MPI_ID==Nodes-1?0:MPI_ID+1);

	if(MPI_ID==0)MPI_Send(&MPI_ID,1,MPI_INT,nextNode,0,MPI_COMM_WORLD);
	int receive_msg;
	MPI_Status recv_state;
	MPI_Recv(&receive_msg,1,MPI_INT,prevNode,0,MPI_COMM_WORLD,&recv_state);

	cout <<"Greeting from Node"<<MPI_ID<<"/"<<Nodes<<" at "<<name<<endl;

	if(MPI_ID!=0)MPI_Send(&MPI_ID,1,MPI_INT,nextNode,0,MPI_COMM_WORLD);
}

/*****INIT OMP*****/
void RunManager::Init_Parallel(){

	//for omp, workspace start at 0; for mpi, workspace start at SWAP_LEN
	int StartX;
	Nodes>1?StartX=SWAP_LEN:StartX=0;

	int processes;

#pragma omp parallel
	{
#pragma omp single
		processes=omp_get_num_threads();
	}

	cout<<"Node "<<MPI_ID<<" has "<<processes<<" threads"<< endl;

#if(CUT==CUTX)
	if(grid->workLength()%processes!=0)cerr<<"Warning: Grid size does not match CPU architecture!"<<endl;
	int divide=(grid->workLength())/processes;
#elif(CUT==CUTY)
	if(grid->gridY()%processes!=0)cerr<<"Warning: Grid size does not match CPU architecture!"<<endl;
	int divide=grid->gridY()/processes;
#elif(CUT==CUTZ)
	if(grid->gridZ()%processes!=0)cerr<<"Warning: Grid size does not match CPU architecture!"<<endl;
	int divide=grid->gridZ()/processes;
#endif

	ranges=new Range[processes]();

#pragma omp parallel
	{
		int i=omp_get_thread_num();
#if(CUT==CUTX)
		ranges[i]=Range(StartX+divide*i,StartX+(i+1)*divide,0,grid->gridY(),0,grid->gridZ());
#elif(CUT==CUTY)
		ranges[i]=Range(StartX,grid->gridX()-StartX,i*divide,(i+1)*divide,0,grid->gridZ());
#elif(CUT==CUTZ)
		ranges[i]=Range(StartX,grid->gridX()-StartX,0,grid->gridY(),i*divide,(i+1)*divide);
#endif
	}

}


