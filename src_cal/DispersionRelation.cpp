/*
 * DispersionRelation.cpp
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#include <Cell.h>
#include <CSPIC.h>
#include <DispersionRelation.h>
#include <Electron.h>
#include <Grid.h>
#include <LePlasma.h>
#include <MPIGrid.h>
#include <PhysicalParameters.h>
#include <RunManager.h>
#include <stdlib.h>
#include <Vector3D.h>
#include <Vertex.h>
#include <cmath>
#include <iostream>

DispersionRelation::DispersionRelation(){

//	You should enable "updateA(range)" and "updateY(range)".
	deltaT=1/(10*LIGHT_SPEED);

	if(RunManager::Nodes>1){
		grid = new MPIGrid(4,1,1);
	}else{
		grid = new Grid(32,1,1,true);
	}

	particle=new Electron();
	particleCount=32*1*1*100;

	thermalVelocity=0.001*LIGHT_SPEED;//0.1*LIGHT_SPEED;

	aVx=aVy=aVz=0;

	//select Engine

	engine=new CSPIC(grid,deltaT);

	launch(REPORT);
}


//void DispersionRelation::report(){
//	if(REPORT){
//		Vertex* v;
//
//		//for_each_Vertex_within(grid,v,grid->workSpace){
//		for_each_Vertex(grid,v){
//
//			outputFile<<v->Y.x<<'\t';
//
//		}end_for_each_Vertex_within
//
//		outputFile<<'\n';
//	}
//}
void DispersionRelation::report(){
	if(REPORT){
		outputFile<<grid->vertex(0,0,0)->Y.x << '\t';
	}
}


DispersionRelation::~DispersionRelation() {
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

void DispersionRelation::distributeParticle(){

	double sigma = thermalVelocity / sqrt(2);

	Range r=grid->World;

	for (int i = 0; i < particleCount; i++){
		Particle *newParticle = particle->clone();

		newParticle->Position.x=Random(r.sx,r.ex);
		newParticle->Position.y=Random(r.sy,r.ey);
		newParticle->Position.z=Random(r.sz,r.ez);

		newParticle->Momentum.x=gaussRand(aVx, sigma);
		newParticle->Momentum.y=gaussRand(aVy, sigma);
		newParticle->Momentum.z=gaussRand(aVz, sigma);

		grid->directAddParticle(newParticle);
	}

}


void DispersionRelation::initP(){
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


void DispersionRelation::initA(){
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


void DispersionRelation::initY(){
	//	Particle *curParticle;
	//	Vector3D r;  //pointing from source point to field point
	//	double lengthR;  //the length of r
	//
	//
	//	Vertex *curVertex;
	//	Vector3D VertexRealPosition;
	//	for_each_Particle_within(grid, curParticle, grid->workSpace){
	//		for_each_Vertex_around(grid, curVertex, curParticle, VertexRealPosition){
	//			r = VertexRealPosition - curParticle->Position;
	//			lengthR = sqrt(Square(r.x) + Square(r.y) + Square(r.z));
	//			curVertex->Y = r / (-Cube(lengthR));
	//		}end_for_each_Vertex_around
	//	}end_for_each_Particle(curParticle)


	for (int i = 0; i < grid->workLength(); ++i) {
		for (int j = 0; j < grid->gridY(); ++j){
			for (int k = 0; k < grid->gridZ(); ++k){
				double worldX=i+grid->workSpace.sx+grid->workLength()*RunManager::MPI_ID;
				grid->vertex(i+grid->workSpace.sx,j,k)->Y=Vector3D(0.000005*cos(worldX/grid->World.rangeX()*2*M_PI),0,0);
			//	grid->vertex(i+grid->workSpace.sx,j,k)->Y=Vector3D(0.000005,0,0);
			}
		}
	}
}
