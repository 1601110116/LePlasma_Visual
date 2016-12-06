/*
 * Distributer.h
 *
 *  Created on: Apr 6, 2016
 *      Author: zlstudio
 */

#ifndef DISTRIBUTER_H_
#define DISTRIBUTER_H_

#include <Particle.h>
#include <Cell.h>
#include <Grid.h>
#include <Electron.h>
#include <MPIGrid.h>

class Cell;

class Distributer {
	public:
		Distributer(int _routes,int _thisroute,Grid* mananged_grid);
		virtual ~Distributer();

		void pack(Particle *p,int dest);

		void assembleAndDispatch();

	private:
		void dispatch(int dest);

		void collect(int dest);

		Grid* grid;

		int routes;
		int thisRoute;

		double** sendBuffer;
		int* sendLength;

		Cell* particleBuffer;
};

#endif /* DISTRIBUTER_H_ */
