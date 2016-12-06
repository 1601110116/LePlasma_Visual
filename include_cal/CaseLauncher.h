/*
 * CaseLauncher.h

 *
 *  Created on: 2016-4-12
 *      Author: lyt
 */

#ifndef INCLUDE_CAL_CASELAUNCHER_H_
#define INCLUDE_CAL_CASELAUNCHER_H_

#include <iostream>
#include <Grid.h>
#include <MPIGrid.h>
#include <Particle.h>
#include <stdlib.h>
#include <Electron.h>
#include <math.h>
#include <PhysicalParameters.h>
#include <CSPIC.h>
#include <Case.h>

using namespace std;

class CaseLauncher {
	public:
		CaseLauncher();
		virtual ~CaseLauncher();

		/* gaussRand: generate a vector according to 3-dimentional gauss distribution function. */
		static double gaussRand(double expectation, double variance);

		/* launchCase: Sets the parameters of the case after initializing the grid. */
		static void launchCase(Grid *grid);
		static void launchCase(Case* _case);
		static void MPILaunchCase(MPIGrid *grid,int totalParticles, int thisID);

		/* distributeParticle: Creats given kind and given amount of particles
		 * and set their initial position and velocity. */
		static void distributeParticle(Grid *grid, Particle *particle, int amount, const Range &r);

		/* setMaxwellVelocity: set the velocity of all particles in grid to Maxwell distribution. */
		static void setMaxwellVelocity(Grid* grid, double thermalVelocity,
						double averageVelocityX, double averageVelocityY, double averageVelocityZ);

		/* initCase: set all A, P, and Y in the grid according to positions and velocities of the particles in grid. */
		static void initA(Grid *grid, const Range &r);
		static void initP(Grid *grid, const Range &r);
		static void initY(Grid *grid, const Range &r);
};



#endif /* INCLUDE_CAL_CASELAUNCHER_H_ */
