/*
 * DispersionRelation.h
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#ifndef INCLUDE_DISPERSIONRELATION_H_
#define INCLUDE_DISPERSIONRELATION_H_

#include <Case.h>

class DispersionRelation:public Case {
	public:
		DispersionRelation();
		virtual ~DispersionRelation();

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
