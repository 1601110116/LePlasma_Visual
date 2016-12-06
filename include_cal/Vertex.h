/*
 * Vertex.h
 *
 *  Created on: Mar 4, 2016
 *      Author: zlstudio
 */

#ifndef VERTEX_H_
#define VERTEX_H_

#include "Vector3D.h"

class Vertex {
	public:
		Vertex();
		Vertex(int x,int y,int z);

		virtual ~Vertex();

		Vector3D A;
		Vector3D Y;

		int x();
		int y();
		int z();

		Vector3D &r();

		string toString();

		int _x,_y,_z;

		Vector3D _r;

};

#endif /* VERTEX_H_ */
