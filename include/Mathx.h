/*
 * Mathx.h
 *
 *  Created on: Apr 23, 2016
 *      Author: zlstudio
 */

#ifndef SRC_MATHX_H_
#define SRC_MATHX_H_

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

struct Arrange{
	double start,end,tick,base;
};

class Mathx {
	public:
		Mathx();
		virtual ~Mathx();

		static Arrange Tick125(double s,double e);

		static double ceilFloat(double n,int d);

		static double floorFloat(double n,int d);

		static string showFloat(double n,int d);

		static string DoubleToString(double d);

};

#endif /* SRC_MATHX_H_ */
