/*
 * Eular.h
 *
 *  Created on: Mar 6, 2016
 *      Author: zlstudio
 */

#ifndef EULAR_H_
#define EULAR_H_

#include "Range.h"
#include "Grid.h"


class Eular {
	public:
		static void nextStep(const Range &,double dt);
		static void nextStep2(const Range &,double dt);
		static void nextStep3(const Range &,double dt);
		static void nextStep4(const Range &,double dt);

		static void updateX(const Range &,double dt);
		static void updateP(const Range &,double dt);
		static void updateY(const Range &,double dt);
		static void updateA(const Range &,double dt);

		static Grid* grid;
};

#endif /* EULAR_H_ */
