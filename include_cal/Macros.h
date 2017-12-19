/*
 * PhysicalParameters.h
 *
 *  Created on: 2016-4-14
 *      Author: lyt
 */

#ifndef INCLUDE_PHYSICALPARAMETERS_H_
#define INCLUDE_PHYSICALPARAMETERS_H_

/* We use Gauss cgs units to manipulate CSPIC. In our simulation we use our DIY units, which we call
 * it Grid units. In Grid units we set
 *         DELTA_X = 2.4355e-2 cm                              which is the size of grid
 *         ME = 9.1094e-28 g                                       which is the mass of a electron
 *         QE = 4.0832e-10 g^0.5cm^1.5s^-1         which is the charge of a eletron
 * to be basic units, so that
 *         g = 1.0978e27 ME
 *         cm = 4.1059e1 DELTA_X
 *         s = 4.1870e6 M_E^0.5DELTA_X^1.5QE^-1
 * so that in Grid units, we have
 *         c = 2.9398e5 M_E^-0.5DELTA_X^-0.5QE
 *         1VOLT = 1.6914e5 DELTA_X^-1QE
 *         1METER = 4.1059e3 DELTA_X
 *         1V/M = 4.1194e1 DELTA_X^-2QE
 *
 * When we input a parameter, we should multiply it with its unit. When a parameter is output, it
 * shoud be devided by its unit.
 *
 */

//#define DELTA_X 1
//#define ME 1
//#define QE 1
//
//#define G 1.0978e27
//#define CM 41.059
//#define SECOND 4.187e6
//
//#define LIGHT_SPEED 2.9398e5
//#define METER 4105.9
//#define VOLT 1.6914e5
//#define KV 1.6914e8
//
#define METER 1.0
#define FOUR_PI 12.56637061



//when DELTA_X = 9.3802e-4 cm

//#define G 1.0978e27
//#define CM 1066.1
//#define SECOND 5.5395e8
//
//#define METER 1.06607e5
//#define LIGHT_SPEED 57676.8
//#define VOLT 6514.27
//#define KV 6.51427e6


//when DELTA_X = 1e-4 cm
//
//#define G 1.0978e27
//#define CM 1e4
//#define SECOND 1.3529e10
//
//#define LIGHT_SPEED 22159.0


//a small case
//#define G 3.13e23
//#define CM 41.059
//#define SECOND 2.1e8
//
//#define LIGHT_SPEED 5866.0



#define Square(a) (a)*(a)
#define Cube(a) (a)*(a)*(a)

//ITER and EAST
#define LIGHT_SPEED 3.2151e1



#endif /* INCLUDE_PHYSICALPARAMETERS_H_ */
