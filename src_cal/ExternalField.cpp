/*
 * ExternalField.cpp
 *
 *  Created on: Mar 13, 2016
 *      Author: zlstudio
 */

#include "ExternalField.h"

ExternalField::ExternalField() {
	// TODO Auto-generated constructor stub

}

ExternalField::~ExternalField() {
	// TODO Auto-generated destructor stub
}

void ExternalField::setExternalField(Grid* grid){

	Vertex *vertex_iterator;

	for_each_Vertex(grid,vertex_iterator){
		if(_i>grid->gridX()/2){
			vertex_iterator->A=Vector3D(-0.1,0,0);
		}else{
			vertex_iterator->A=Vector3D(0.1,0,0);
		}
	}end_for_each_Vertex
}
