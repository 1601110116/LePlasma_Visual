/*
 * UniformE.cpp
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#include <Cell.h>
#include <CSPIC.h>
#include <UniformE.h>
#include <Electron.h>
#include "Grid.h"
#include <LePlasma.h>
#include <MPIGrid.h>
#include "PhysicalParameters.h"
#include <RunManager.h>
#include <stdlib.h>
#include <Vector3D.h>
#include <Vertex.h>
#include <cmath>
#include <iostream>

UniformE::UniformE(){

	deltaT=1/(10*LIGHT_SPEED);

	if(RunManager::Nodes>1){
		grid = new MPIGrid(120,1,1);
	}else{
		grid = new Grid(120,1,1,true);
	}

	particle=new Electron();
	particleCount=1;

	thermalVelocity=0;//0.1*LIGHT_SPEED;

	aVx=aVy=aVz=0;

	//select Engine

	engine=new CSPIC(grid,deltaT);

	launch(REPORT);
}


void UniformE::report(){
	if(REPORT){
		Particle *curParticle;
		for_each_Particle(grid,curParticle){
			outputFile << curParticle->Position.x << '\t';
		}end_for_each_Particle(curParticle)
	}
}

UniformE::~UniformE() {
	// TODO Auto-generated destructor stub
}


inline double Random(double min, double max){
	return (rand()*(max-min))/RAND_MAX+min;
}

static double gaussRand(double expectation, double sigma){
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

void UniformE::distributeParticle(){

	double sigma = thermalVelocity / sqrt(2);

	Range r=grid->World;

	for (int i = 0; i < particleCount; i++){
		Particle *newParticle = particle->clone();

		newParticle->Position.x=0;
		newParticle->Position.y=0.5;
		newParticle->Position.z=0.5;

		newParticle->Momentum.x=gaussRand(aVx, sigma);
		newParticle->Momentum.y=gaussRand(aVy, sigma);
		newParticle->Momentum.z=gaussRand(aVz, sigma);

		grid->directAddParticle(newParticle);
	}

}


void UniformE::initP(){
	Particle *curParticle;
	Vertex *curVertex;

	//Vector3D v1, v2,
	Vector3D VertexRealPosition;
	for_each_Particle_within(grid, curParticle,grid->workSpace){
		//currently Momentum is Velocity
		for_each_Vertex_around(grid, curVertex, curParticle,VertexRealPosition){
			curParticle->Momentum -= curVertex->A * _W(curParticle->Position - VertexRealPosition);
		}end_for_each_Vertex_around

	}end_for_each_Particle(curParticle)
}


void UniformE::initA(){
	//	Particle *curParticle;
	//	double rc2;  //the distance between the source point and the field point multimplied by c square
	//
	//	Vertex *curVertex;
	//	Vector3D VertexRealPosition;
	//	Vector3D r;
	//	for_each_Particle_within(grid, curParticle,grid->workSpace){
	//		//currently Momentum is Velocity
	//		for_each_Vertex_around(grid, curVertex, curParticle, VertexRealPosition){
	//			r = VertexRealPosition - curParticle->Position;
	//			rc2 = sqrt(Square(r.x) + Square(r.y) + Square(r.z)) * Square(LIGHT_SPEED);
	//			curVertex->A += (curParticle->Momentum / rc2);
	//		}end_for_each_Vertex_around
	//
	//	}end_for_each_Particle(curParticle)
}


void UniformE::initY(){
	for (int i = 0; i < grid->gridX(); ++i) {
		for (int j = 0; j < grid->gridY(); ++j){
			for (int k = 0; k < grid->gridZ(); ++k){
				grid->vertex(i,j,k)->Y=\
						Vector3D(240.0/(Square(100.0/LIGHT_SPEED)*4.0*M_PI),0,0);
			}
		}
	}
}
