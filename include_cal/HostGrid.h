/*
 * HostGrid.h
 *
 *  Created on: Apr 13, 2016
 *      Author: zlstudio
 */

#ifndef HOSTGRID_H_
#define HOSTGRID_H_


#if REPORT
#include <inet.h>
#endif


#include "Grid.h"

class HostGrid : public Grid{
	public:
		HostGrid(int x,int y,int z,int n);
		virtual ~HostGrid();

		void askReport();

	private:

		int nodes;

		int workLength;
#if REPORT
		int reportSocket;
#endif

		Range* space;

		bool receiveParticlesFrom(int src);

		bool receiveVertexesFrom(int src);
};

#endif /* HOSTGRID_H_ */
