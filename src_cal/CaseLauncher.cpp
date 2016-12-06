/*
 * CaseLauncher.cpp

 *
 *  Created on: 2016-4-12
 *      Author: lyt
 */

#include <cmath>
#include <cstdlib>
#include "MPIGrid.h"
#include "PhysicalParameters.h"
#include "Grid.h"
#include "Case.h"
#include "CaseLauncher.h"
#include "RunManager.h"
#include "LePlasma.h"
#include "Electron.h"
#include "Cell.h"
#include "Vertex.h"
#include "CSPIC.h"

inline double Random(double min, double max){
	return (rand()*(max-min))/RAND_MAX+min;
}

void CaseLauncher::launchCase(Case* _case){
	srand(0);

	if(RunManager::MPI_ID==0){
		_case->distributeParticle();
	}

	_case->grid->refreshParticleLocation();

	cout<<"Node "<<RunManager::MPI_ID<<" Distribute Complete. Particles: "<<_case->grid->particles()<<endl;

	_case->initA();
	_case->initP();
	_case->initY();

	if(RunManager::Nodes>1){
		_case->grid->refreshParticleLocation();
	}
}

void CaseLauncher::launchCase(Grid *grid){
	srand(0);

	distributeParticle(grid, new Electron(), 120*100, grid->workSpace);

	setMaxwellVelocity(grid,0.031267*LIGHT_SPEED,0,0,0);

	grid->refreshParticleLocation();

	//setMaxwellVelocity(grid,0.02513*LIGHT_SPEED,0,0,0);

	//setMaxwellVelocity(grid,0,1000,0,0);

	//	initA(grid, grid->workSpace);
	initP(grid,grid->workSpace);
	//	initY(grid, grid->workSpace);
}

void CaseLauncher::MPILaunchCase(MPIGrid *grid,int totalParticles,int thisID){
	srand(0);

	if(thisID==0){
		distributeParticle(grid, new Electron(), totalParticles,grid->World);
		setMaxwellVelocity(grid,0.031267*LIGHT_SPEED,0,0,0);
	}
	grid->refreshParticleLocation();

	//	initA(grid, grid->workSpace);
	initP(grid,grid->workSpace);
	//	initY(grid, grid->workSpace);

	grid->refreshParticleLocation();
}

CaseLauncher::CaseLauncher(){

}
CaseLauncher::~CaseLauncher(){

}

double CaseLauncher::gaussRand(double expectation, double sigma){
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if (phase == 0) {
		do{
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		}  while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	X = X * sigma + expectation;
	return X;
}

#define Reverse true

void  CaseLauncher::setMaxwellVelocity(Grid* grid, double thermalVelocity, double averageVelocityX, double averageVelocityY, double averageVelocityZ){

	double sigma = thermalVelocity / sqrt(2);

	Particle* p;
#if Reverse
	Vector3D tmp[grid->particles()];

	for (int i = 0; i < grid->particles(); ++i) {
		tmp[i].x=gaussRand(averageVelocityX, sigma);
		tmp[i].y=gaussRand(averageVelocityY, sigma);
		tmp[i].z=gaussRand(averageVelocityZ, sigma);
	}

	int n=0;

	for_each_Particle_within(grid,p,grid->workSpace){

		p->Momentum.x=tmp[grid->particles()-n-1].x;
		p->Momentum.y=tmp[grid->particles()-n-1].y;
		p->Momentum.z=tmp[grid->particles()-n-1].z;

		n++;

	}end_for_each_Particle(p)

#else
	for_each_Particle_within(grid,p,grid->workSpace){

		p->Momentum.x=gaussRand(averageVelocityX, sigma);
		p->Momentum.y=gaussRand(averageVelocityY, sigma);
		p->Momentum.z=gaussRand(averageVelocityZ, sigma);

	}end_for_each_Particle(p)


#endif

}

void CaseLauncher::distributeParticle(Grid *grid, Particle *particle, int amount ,const Range &r){
	for (int i = 0; i < amount; i++){
		Particle *newParticle = particle->clone();
		newParticle->Position = Vector3D(
				Random(r.sx,r.ex),
				Random(r.sy,r.ey),
				Random(r.sz,r.ez)
		);
		grid->directAddParticle(newParticle);
	}
}

void CaseLauncher::initA(Grid *grid, const Range &range){
	Particle *curParticle;
	double rc2;  //the distance between the source point and the field point multimplied by c square

	Vertex *curVertex;
	Vector3D VertexRealPosition;
	Vector3D r;
	for_each_Particle_within(grid, curParticle,range){
		//currently Momentum is Velocity
		for_each_Vertex_around(grid, curVertex, curParticle, VertexRealPosition){
			r = VertexRealPosition - curParticle->Position;
			rc2 = sqrt(Square(r.x) + Square(r.y) + Square(r.z)) * Square(LIGHT_SPEED);
			curVertex->A += (curParticle->Momentum / rc2);
		}end_for_each_Vertex_around

	}end_for_each_Particle(curParticle)


	//	Vertex *curVertex;
	//	for_each_Vertex_within(grid, curVertex, r){
	//		curVertex->A = Vector3D(-curVertex->y()/2.0, curVertex->x()/2.0, 0)*(-LIGHT_SPEED/50);
	//	}end_for_each_Vertex_within
}

void CaseLauncher::initP(Grid *grid, const Range &r){
	Particle *curParticle;
	Vertex *curVertex;

	//Vector3D v1, v2,
	Vector3D VertexRealPosition;
	for_each_Particle_within(grid, curParticle,r){
		//currently Momentum is Velocity
		for_each_Vertex_around(grid, curVertex, curParticle,VertexRealPosition){
			curParticle->Momentum -= curVertex->A * _W(curParticle->Position - VertexRealPosition);
		}end_for_each_Vertex_around

	}end_for_each_Particle(curParticle)
}

void CaseLauncher::initY(Grid *grid,const Range &range){
	Particle *curParticle;
	Vector3D r;  //pointing from source point to field point
	double lengthR;  //the length of r


	Vertex *curVertex;
	Vector3D VertexRealPosition;
	for_each_Particle_within(grid, curParticle, range){
		for_each_Vertex_around(grid, curVertex, curParticle, VertexRealPosition){
			r = VertexRealPosition - curParticle->Position;
			lengthR = sqrt(Square(r.x) + Square(r.y) + Square(r.z));
			curVertex->Y = r / (-Cube(lengthR));
		}end_for_each_Vertex_around
	}end_for_each_Particle(curParticle)



	//	int i,j,k;
	//	for(i = 0; i < grid->gridX(); i++){
	//		for (j = 0; j < grid->gridZ(); j++){
	//			grid->vertex(i,0,j)->Y.x = 10;
	//		}
	//	}



	//	double K = M_PI/2;
	//	Vertex *curVertex;
	//	for_each_Vertex_within(grid, curVertex, range){
	//		curVertex->Y.x = -10*cos(K*curVertex->r().x);
	//	}end_for_each_Vertex_within
}
