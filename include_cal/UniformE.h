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

class UniformE:public Case {
	public:
		UniformE();
		virtual ~UniformE();

		void distributeParticle();
		void initP();
		void initA();
		void initY();

		void report();

		Particle* particle;
		int particleCount;

		double thermalVelocity;
		double aVx,aVy,aVz;
};

#endif /* INCLUDE_DISPERSIONRELATION_H_ */
