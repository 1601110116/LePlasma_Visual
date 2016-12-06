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

class UniformB:public Case {
	public:
		UniformB();
		virtual ~UniformB();

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
