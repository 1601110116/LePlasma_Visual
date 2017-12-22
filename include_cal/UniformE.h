/*
 * UniformE.h
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#ifndef INCLUDE_UNIFORME_H_
#define INCLUDE_UNIFORME_H_

#include "Particle.h"
#include "Case.h"
#include <map>

class UniformE:public Case {
public:
	UniformE();
	virtual ~UniformE();

	void distributeParticle() override;
	void initP() override;
	void initA() override;
	void initY() override;

	void report() override;

	Particle* particle;
	int particleCount;
	map<string, double> units;
	double lightSpeed;

	double thermalVelocity;
	double aVx,aVy,aVz;
};

#endif /* INCLUDE_DISPERSIONRELATION_H_ */
