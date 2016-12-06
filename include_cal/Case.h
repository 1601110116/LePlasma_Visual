/*
 * Case.h
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#ifndef SRC_CASE_H_
#define SRC_CASE_H_

#include <iostream>
#include <fstream>


class Grid;
class Engine;
class Range;

using namespace std;


class Case {
	public:
		Case();
		virtual ~Case();

		virtual void distributeParticle();
		virtual void initP();
		virtual void initA();
		virtual void initY();
		virtual void report();

		virtual void launch(bool report);

		Grid* grid;

		double deltaT;

		Engine* engine;

		ofstream outputFile;

};

#endif /* SRC_CASE_H_ */
