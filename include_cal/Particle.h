/*
 * Particle.h
 *
 *  Created on: Mar 3, 2016
 *      Author: zlstudio
 */

#ifndef PARTICAL_H_
#define PARTICAL_H_

#include "Vector3D.h"
#include <Cell.h>
#include <PhysicalParameters.h>

class Cell;
class Vector3D;

class Particle {
	public:
		Particle();
		virtual ~Particle();

		double charge;
		double mass;
		const char* name;

		Vector3D Position;
		Vector3D Momentum;

		Cell* cell;

		Particle* nextParticle;
		Particle* prevParticle;

		virtual Particle* clone();

#if USE_CACHE
//A 4*4*4 matrix to store cache

		double W_cache[4][4][4];
		Vector3D GW_cache[4][4][4];

#endif

	protected:
		Particle(double,double,const char*);

};

#endif /* PARTICAL_H_ */
