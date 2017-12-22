/*
 * UniformB.cpp
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#include <Cell.h>
#include <UniformB.h>
#include <Electron.h>
#include <Grid.h>
#include <LePlasma.h>
#include <MPIGrid.h>
#include <Macros.h>
#include <RunManager.h>
#include <stdlib.h>
#include <Vector3D.h>
#include <Vertex.h>
#include <cmath>
#include <iostream>
#include "EngineForSingleParticle.h"
#include "CSPIC.h"



inline double _W1(double x){
	if (x > 2)
		return 0.0;
	else if (x > 1)
		return x * (Cube(x) * (x * (x * (x * (15.0/1024 * x - 15.0/128) + 49.0/128) - 21.0/32) + 35.0/64) - 1.0) + 1.0;
	else if (x > 0)
		return Square(x) * (Square(x) * (x * (x * (x * (-15.0/1024 * x - 15.0/128) + 7.0/16) - 21.0/32) + 175.0/256) - 105.0/128) + 337.0/512;
	else if (x > -1)
		return x * x * (x * x * (x * (x * (x * (-15.0/1024 * x + 15.0/128) + 7.0/16) + 21.0/32) + 175.0/256) - 105.0/128) + 337.0/512;
	else if (x > -2)
		return x * (Cube(x) * (x * (x * (x * (15.0/1024 * x + 15.0/128) + 49.0/128) + 21.0/32) + 35.0/64) + 1.0) + 1.0;
	else
		return 0.0;
}

inline double _W(const Vector3D &r){
	return _W1(r.x)*_W1(r.y)*_W1(r.z);
}

void UniformB::calcUnits() {
    // cell width in cm
    double dx = 5.8594e-2;
    // super particle weight in g
    double mp = 1.8324e-19;
    // super particle charge in g^1/2cm^3/2s^-1
    double qp = 9.6616e-2;

    units["cm"] = 1.0/dx;
    units["gram"] = 1.0/mp;
    units["second"] = qp*pow(units["gram"],0.5)*pow(units["cm"],1.5);
    units["lightSpeed"] = 2.9979e10*units["cm"]/units["second"];
    units["gauss"] = pow(units["gram"],0.5) * pow(units["cm"],-0.5) * pow(units["second"],-1);
    units["omegaCi"] = 9.58e3*1.0e-2/units["second"];

    lightSpeed = units["lightSpeed"];

}

UniformB::UniformB(){
    calcUnits();

/*	You should disable "updateA(range)" and "updateY(range)" first */
	deltaT=5.0*M_PI/(lightSpeed);

	if(RunManager::Nodes>1){
		grid = new MPIGrid(20,20,1);
	}else{
		grid = new Grid(20,20,1,true);
	}

	particle=new Electron();
	particleCount=1;

	thermalVelocity=0;//0.1*LIGHT_SPEED;

	aVx=aVy=aVz=0;
	aVy=-0.1*lightSpeed;

	//select Engine

	engine=new CSPIC(grid,deltaT,units);

	launch(REPORT);
}

void UniformB::launch(bool report){

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

void UniformB::report(){
	if(REPORT){
		Particle *curParticle;
		for_each_Particle(grid,curParticle){
			outFileX << curParticle->Position.x << '\t';
			outFileY << curParticle->Position.y << '\t';
		}end_for_each_Particle(curParticle)
	}
}


UniformB::~UniformB() {
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

void UniformB::distributeParticle(){

	double sigma = thermalVelocity / sqrt(2);

	Range r=grid->World;

	for (int i = 0; i < particleCount; i++){
		Particle *newParticle = particle->clone();

		newParticle->Position.x=5;
		newParticle->Position.y=10;
		newParticle->Position.z=0.5;

		newParticle->Momentum.x=gaussRand(aVx, sigma);
		newParticle->Momentum.y=gaussRand(aVy, sigma);
		newParticle->Momentum.z=gaussRand(aVz, sigma);

		grid->directAddParticle(newParticle);
	}

}


void UniformB::initP(){
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


void UniformB::initA(){
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
		curVertex->A=Vector3D(-curVertex->y()/2.0,curVertex->x()/2.0,0)*(lightSpeed/50.0);
	}end_for_each_Vertex_within
}


void UniformB::initY(){
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
