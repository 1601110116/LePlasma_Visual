/*
 * Electron.h
 *
 *  Created on: Mar 4, 2016
 *      Author: zlstudio
 */

#ifndef ELECTRON_H_
#define ELECTRON_H_

#include "Particle.h"

class Particle;

class Electron:public Particle {
	public:
		Electron();
		virtual ~Electron();

		Electron* clone();

};

#endif /* ELECTRON_H_ */
