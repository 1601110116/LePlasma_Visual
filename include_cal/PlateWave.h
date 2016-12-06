/*
 * PlateWave.h
 *
 *  Created on: 2016-5-9
 *      Author: zlstudio
 */

#ifndef SRC_CAL_PLATEWAVE_H_
#define SRC_CAL_PLATEWAVE_H_

class Grid;
class Range;

class PlateWave {
	public:
		PlateWave();
		virtual ~PlateWave();

		static void launchCase(Grid *grid);

		static void initY(Grid *grid, const Range &r);
};

#endif /* SRC_CAL_PLATEWAVE_H_ */
