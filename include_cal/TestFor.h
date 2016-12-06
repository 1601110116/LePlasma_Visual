/*
 * DispersionRelation.h
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#ifndef INCLUDE_TESTFOR_H_
#define INCLUDE_TESTFOR_H_

#include "Case.h"
#include "Particle.h"

class TestFor:public Case {
	public:
		TestFor();
		virtual ~TestFor();

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
