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
#include <Macros.h>

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
	Vector3D A;
//	Vector3D B;  // the local magnetic field in grid units
	Vector3D X;  // the position of the particle in cartesian coordinate in grid units

	Cell* cell;

	Particle* nextParticle;
	Particle* prevParticle;

	virtual Particle* clone();

//A 4*4*4 matrix to store cache

	double W_cache[4][4][4];
	Vector3D GW_cache[4][4][4];


protected:
	Particle(double,double,const char*);

};

#endif /* PARTICAL_H_ */
