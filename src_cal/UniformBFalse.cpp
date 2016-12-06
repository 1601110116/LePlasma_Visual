/*
 * UniformBFalse.cpp
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#include <Cell.h>
#include <CSPIC.h>
#include <UniformBFalse.h>
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

UniformBFalse::UniformBFalse(){

/*	You should change "updateP(range)" into "updatePFalse(range)" first,
 *     and then disable "updateA(range)" and "updateY(range)".
 */
	deltaT=M_PI/(10*LIGHT_SPEED);

	if(RunManager::Nodes>1){
		grid = new MPIGrid(130,130,1);
	}else{
		grid = new Grid(130,130,1,true);
	}

	particle=new Electron();
	particleCount=1;

	thermalVelocity=0;//0.1*LIGHT_SPEED;

	aVx=aVy=aVz=0;
	aVy=-0.1*LIGHT_SPEED;

	//select Engine

	engine=new CSPIC(grid,deltaT);

	launch(REPORT);
}

void UniformBFalse::launch(bool report){

	srand(0);

	if(RunManager::MPI_ID==0){
		distributeParticle();
	}

	grid->refreshParticleLocation();

	cout<<"Node "<<RunManager::MPI_ID<<" Distribute Complete. Particles: "<<grid->particles()<<endl;

	initA();
	initP();
	initY();

	if(RunManager::Nodes>1){
		grid->refreshParticleLocation();
	}

	if(report){
		time_t t = time(0);
		char tmp[64];
		strftime(tmp, sizeof(tmp), "%Y_%m_%d_%H_%M",localtime(&t));
		ostringstream name1;
		ostringstream name2;
		name1<<"Report"<<tmp<<"PositionX" << ".output";
		name2 << "Report" << tmp << "PositionY" << ".output";
		outFileX.open(name1.str().c_str());
		outFileY.open(name2.str().c_str());
	}
}

void UniformBFalse::report(){
	if(REPORT){
		Particle *curParticle;
		for_each_Particle(grid,curParticle){
			outFileX << curParticle->Position.x << '\t';
			outFileY << curParticle->Position.y << '\t';
		}end_for_each_Particle(curParticle)
	}
}


UniformBFalse::~UniformBFalse() {
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

void UniformBFalse::distributeParticle(){

	double sigma = thermalVelocity / sqrt(2);

	Range r=grid->World;

	for (int i = 0; i < particleCount; i++){
		Particle *newParticle = particle->clone();

		newParticle->Position.x=60;
		newParticle->Position.y=65;
		newParticle->Position.z=0.5;

		newParticle->Momentum.x=gaussRand(aVx, sigma);
		newParticle->Momentum.y=gaussRand(aVy, sigma);
		newParticle->Momentum.z=gaussRand(aVz, sigma);

		grid->directAddParticle(newParticle);
	}

}


void UniformBFalse::initP(){
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


void UniformBFalse::initA(){
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
	Vertex *curVertex;
	for_each_Vertex_within(grid,curVertex,grid->workSpace){
		curVertex->A=Vector3D(-curVertex->y()/2.0,curVertex->x()/2.0,0)*(LIGHT_SPEED/50.0);
	}end_for_each_Vertex_within
}


void UniformBFalse::initY(){
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


}
