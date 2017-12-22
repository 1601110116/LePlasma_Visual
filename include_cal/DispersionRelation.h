/*
 * DispersionRelation.h
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#ifndef INCLUDE_DISPERSIONRELATION_H_
#define INCLUDE_DISPERSIONRELATION_H_

#include <Case.h>
#include "Particle.h"
#include <map>

class DispersionRelation:public Case {
public:
	DispersionRelation();
	virtual ~DispersionRelation();

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
