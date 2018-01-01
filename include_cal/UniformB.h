/*
 * UniformB.h
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#ifndef INCLUDE_UNIFORMB_H_
#define INCLUDE_UNIFORMB_H_

#include "Case.h"
#include "Particle.h"
#include <map>

class UniformB:public Case {
public:
	UniformB();
	virtual ~UniformB();

	void distributeParticle() override;
	void initP() override;
	void initA() override;
	void initY() override;
	void launch(bool report) override;
    void calcUnits() override ;

	void report() override;

	Particle* particle;
	int particleCount;
	double lightSpeed;

	double thermalVelocity;
	double aVx,aVy,aVz;
	ofstream outFileX;
	ofstream outFileY;
};

#endif /* INCLUDE_DISPERSIONRELATION_H_ */
