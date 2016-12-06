/*
 * Particle.cpp
 *
 *  Created on: Mar 5, 2016
 *      Author: zlstudio
 */

#include <Particle.h>

Particle::Particle():
charge(0),mass(0),name(""),
Position(),Momentum(),
cell(NULL),
nextParticle(NULL),prevParticle(NULL)
{};

Particle::Particle(double c,double m,const char* n):
charge(c),mass(m),name(n),
Position(),Momentum(),
cell(NULL),
nextParticle(NULL),prevParticle(NULL)
{}

Particle::~Particle(){};

Particle* Particle::clone(){
	return new Particle();
};



