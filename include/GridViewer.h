/*
 * GridViewer.h
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#ifndef INCLUDE_GRIDVIEWER_H_
#define INCLUDE_GRIDVIEWER_H_

#include <Cell.h>
#include <GL/freeglut_std.h>
#include <GL/gl.h>
#include <Sprite.h>
#include <Grid.h>

class Grid;



class GridViewer:public Sprite {
	public:
		GridViewer(Grid* grid);
		virtual ~GridViewer();

		void refreshGrid();

		double cellSize;

	private:
		Grid* grid;
		void drawParticles();
		void drawGrid();
};

#endif /* INCLUDE_GRIDVIEWER_H_ */
