/*
 * GridViewer.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#include <GridViewer.h>

GridViewer::GridViewer(Grid* grid):grid(grid),cellSize(20) {
	// TODO Auto-generated constructor stub

}

GridViewer::~GridViewer() {
	// TODO Auto-generated destructor stub
}

void GridViewer::refreshGrid(){

	graphics.clear();
	drawGrid();
	drawParticles();
}

void GridViewer::drawParticles(){
	graphics.pointStyle(0,0xff0000,1);

	Particle* particle_iterator;
	Vector3D pos;

	int n=0;
	for_each_Particle(grid,particle_iterator){
		pos=particle_iterator->Position*cellSize;

		graphics.drawPoint(
				pos.x-grid->gridX()*cellSize/2,
				pos.y-grid->gridY()*cellSize/2,
				pos.z-grid->gridZ()*cellSize/2);
		n++;

	}end_for_each_Particle(particle_iterator);

}
void GridViewer::drawGrid(){

	Cell* cell_iterator;

	graphics.lineStyle(0,0x202020,1);

	for_each_Cell(grid,cell_iterator){
		graphics.drawCube(
				cell_iterator->x()*cellSize-grid->gridX()*cellSize/2,
				cell_iterator->y()*cellSize-grid->gridY()*cellSize/2,
				cell_iterator->z()*cellSize-grid->gridZ()*cellSize/2,
				cellSize);
	}end_for_each_Cell


	graphics.lineStyle(0,0x204020,1);

	graphics.drawBox(-grid->gridX()*cellSize/2,
			-grid->gridY()*cellSize/2,
			-grid->gridZ()*cellSize/2,
			grid->gridX()*cellSize,
			grid->gridY()*cellSize,
			grid->gridZ()*cellSize);
}
