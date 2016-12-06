/*
 * Electron.cpp
 *
 *  Created on: Mar 4, 2016
 *      Author: zlstudio
 */

#include "Electron.h"

Electron::Electron():Particle(-1.0,1.0,"e-")
{
}

Electron::~Electron() {
	// TODO Auto-generated destructor stub
}

Electron* Electron::clone(){
	return new Electron();
}
