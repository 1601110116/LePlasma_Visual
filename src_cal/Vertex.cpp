/*
 * Vertex.cpp
 *
 *  Created on: Mar 4, 2016
 *      Author: zlstudio
 */

#include <Vertex.h>

Vertex::Vertex():_x(0),_y(0),_z(0),A(),Y(),_r(){
}

Vertex::Vertex(int x,int y,int z):_x(x),_y(y),_z(z),A(),Y(),_r(x,y,z){
}


Vertex::~Vertex() {
	// TODO Auto-generated destructor stub
}

string Vertex::toString(){
	ostringstream tmpchar;

	tmpchar<<"Vertex at ("<<_x<<","<<_y<<","<<_z<<")";

	return tmpchar.str();
}

int Vertex::x(){
	return _x;
}

int Vertex::y(){
	return _y;
}

int Vertex::z(){
	return _z;
}

Vector3D& Vertex::r(){
	return _r;
}
