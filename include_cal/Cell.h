/*
 * Cell.h
 *
 *  Created on: Mar 3, 2016
 *      Author: zlstudio
 */

#ifndef CELL_H_
#define CELL_H_

#include <Particle.h>
#include <iostream>
#include <cstddef>

#include <string>
#include <sstream>

#define for_each_Particle_in_Cell(cell,iterator) \
		iterator=cell->firstParticle();\
		while(iterator!=NULL){

#define end_for_each_Particle_in_Cell(iterator) iterator=iterator->nextParticle;}

#define end_for_each_Particle_in_Cell_raw }

using namespace std;

class Particle;

class Cell{
	public:
		Cell();
		Cell(int x,int y,int z);

		virtual ~Cell();

		int x();
		int y();
		int z();

		bool pushParticle(Particle* p);
		bool withdrawParticle(Particle* p);

		//to make use of current particle list when sync.
		int overwriteParticle(Particle* p);
		void restore();
		int cut();

		void destroyAllParticlesInside();

		int length();

		Particle* firstParticle();

		string toString();

		string treeView();

		int _x,_y,_z;

	private:

		Particle* head;
		Particle* tail;
		int _length;

		//to make use of current particle list when sync.
		//retore pointer points to the existed particle to be overwrite

		Particle* restorePointer;
		int restoreLength;
};

#endif /* CELL_H_ */
