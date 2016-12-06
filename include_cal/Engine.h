/*
 * Engine.h
 *
 *  Created on: 2016-5-19
 *      Author: zlstudio
 */

#ifndef SRC_ENGINE_H_
#define SRC_ENGINE_H_

class Range;
class Grid;

class Engine {
	public:
		Engine(Grid*,double);
		virtual ~Engine();

		virtual void update(const Range&);

		Grid* grid;

		double deltaT;

};

#endif /* SRC_ENGINE_H_ */
