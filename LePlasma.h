/*
 * LePlasma_Visual.h
 *
 *  Created on: Apr 16, 2016
 *      Author: zlstudio
 */

#ifndef LEPLASMA_H_
#define LEPLASMA_H_

#define RECORD true

#define REPORT false

#define CUTX 0
#define CUTY 1
#define CUTZ 2
#define CUT CUTX

#define MPI_PARALLEL false

#define DEGUB false
#define OPTIMIZE true
#define OPTIMIZE_1D_MPI false

// if width >= 2 and height=length=1, you should switch on
#define OPTIMIZE_1D false

// if width, height, length >=2, you should switch on
#define OPTIMIZE_3D false

#define USE_CACHE true

#endif /* LEPLASMA_H_ */
