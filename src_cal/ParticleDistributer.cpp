/*
 * ParticleDistributer.cpp
 *
 *  Created on: Mar 5, 2016
 *      Author: zlstudio
 */

#include <ParticleDistributer.h>
#include <math.h>
#include <stdlib.h>

inline double Random(double min, double max){
	return (rand()*(max-min))/RAND_MAX+min;
}

void ParticleDistributer::Distributer_Random(Grid* grid, Particle* p, int amount){

	srand(1);

	for (int i = 0; i < amount; ++i) {

		Particle* particle=p->clone();

//		particle->Position=Vector3D(3.1,5,5);

		particle->Position=Vector3D(
						Random(0,grid->gridX()),
						Random(0,grid->gridY()),
						Random(0,grid->gridZ())
				);

//		particle->Momentum=Vector3D(
//				Random(-10,10),
//				Random(-10,10),
//				Random(-10,10)
//		);

		particle->Momentum=Vector3D(0,-10,0);
		//cout<<particle->Position->toString();
		grid->addParticle(particle);

	}
}

void ParticleDistributer::Distributer_Linear(Grid* grid, Particle* p, int amount){

	srand(1);

	for (int i = 0; i < amount; ++i) {

		Particle* particle=p->clone();

//		particle->Position=Vector3D(3.1,5,5);

		particle->Position=Vector3D(
						10,//double(grid->gridX()*i)/amount,
						grid->gridY()/2,
						grid->gridZ()/2
				);


		particle->Momentum=Vector3D(0,0,0);
		grid->addParticle(particle);

	}
}

void ParticleDistributer::MPIDistributer_Random(MPIGrid* grid, Particle* p, int amount, int node){
	srand(node);

	for (int i = 0; i < amount; ++i) {

		Particle* particle=p->clone();

//		particle->Position=Vector3D(3.1,5,5);

		particle->Position=Vector3D(
						Random(2,2+grid->_workLength),
						Random(0,grid->gridY()),
						Random(0,grid->gridZ())
				);

//		particle->Momentum=Vector3D(
//				Random(-10,10),
//				Random(-10,10),
//				Random(-10,10)
//		);

		//particle->Momentum=Vector3D(1,0,0);
		//cout<<particle->Position->toString();
		grid->addParticle(particle);

	}
}


