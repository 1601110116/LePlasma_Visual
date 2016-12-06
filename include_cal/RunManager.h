/*
 * RunManager.h
 *
 *  Created on: Mar 6, 2016
 *      Author: zlstudio
 */

#ifndef RUNMANAGER_H_
#define RUNMANAGER_H_

#include <string>

class Case;
class MPIGrid;
class Range;


class RunManager {
	public:

		RunManager(Case* _case);

		RunManager(MPIGrid*,void (*Method)(const Range &,double),double);

		static Grid* grid;

		static Case* _Case;

		static double deltaT;

		static int MPI_ID,Nodes;

		static void (*Method)(const Range &,double);

		virtual ~RunManager();

		static void stepNext();

		static void MPI_Init_Prarllel();

	private:

		static Range* ranges;

		static void Init_Parallel();

};


#endif /* MPIRUNMANAGER_H_ */
