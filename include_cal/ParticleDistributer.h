/*
 * ParticleDistributer.h
 *
 *  Created on: Mar 5, 2016
 *      Author: zlstudio
 */

#ifndef PARTICALDISTRIBUTER_H_
#define PARTICALDISTRIBUTER_H_

#include "Grid.h"
#include "MPIGrid.h"

using namespace std;

class ParticleDistributer {
	public:
		 static void Distributer_Random(Grid* grid, Particle* p, int amount);
		 static void Distributer_Linear(Grid* grid, Particle* p, int amount);

		 static void MPIDistributer_Random(MPIGrid* grid, Particle* p, int amount, int node);
};

#endif /* PARTICALDISTRIBUTER_H_ */
