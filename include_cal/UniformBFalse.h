/*
 * UniformBFalse.h
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#ifndef INCLUDE_UNIFORMBFALSE_H_
#define INCLUDE_UNIFORMBFALSE_H_

#include "Case.h"
#include "Particle.h"

class UniformBFalse:public Case {
	public:
		UniformBFalse();
		virtual ~UniformBFalse();

		void distributeParticle();
		void initP();
		void initA();
		void initY();
		void launch(bool report);

		void report();

		Particle* particle;
		int particleCount;

		double thermalVelocity;
		double aVx,aVy,aVz;
		ofstream outFileX;
		ofstream outFileY;
};

#endif /* INCLUDE_DISPERSIONRELATION_H_ */
