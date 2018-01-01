/*
 * Engine.h
 *
 *  Created on: 2016-5-19
 *      Author: zlstudio
 */

#ifndef SRC_ENGINE_H_
#define SRC_ENGINE_H_

#include <map>
#include <string>

using namespace std;
class Range;
class Grid;

class Engine {
public:
	Engine(Grid* grid,double dt,map<string,double>& units);
	virtual ~Engine();

	virtual void update(const Range&);

	Grid* grid;

	double deltaT;
	map<string, double> &units;


};

#endif /* SRC_ENGINE_H_ */
