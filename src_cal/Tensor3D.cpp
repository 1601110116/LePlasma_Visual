/*
 * Tensor3D.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: zlstudio
 */

#include "Tensor3D.h"

Tensor3D::Tensor3D() {
}

Tensor3D::Tensor3D(Vector3D x, Vector3D y, Vector3D z){
	this->x=x;
	this->y=y;
	this->z=z;
}

Tensor3D::~Tensor3D() {
	// TODO Auto-generated destructor stub
}

Tensor3D &Tensor3D::operator =(const Tensor3D &a){
	x=a.x;
	y=a.y;
	z=a.z;
	return *this;
}

Tensor3D &Tensor3D::operator *=(double a){
	x*=a;
	y*=a;
	z*=a;
	return *this;
}

Tensor3D &Tensor3D::operator /=(double a){
	x/=a;
	y/=a;
	z/=a;
	return *this;
}

Tensor3D Tensor3D::operator /(double a)const{
	return Tensor3D(x/a,y/a,z/a);
}

void Tensor3D::restore(){
	x.restore();
	y.restore();
	z.restore();
}
string Tensor3D::toString(){

	ostringstream tmpchar;

	tmpchar<<"|"<<x.x<<","<<x.y<<","<<x.z<<"|"<<endl;
	tmpchar<<"|"<<y.x<<","<<y.y<<","<<y.z<<"|"<<endl;
	tmpchar<<"|"<<z.x<<","<<z.y<<","<<z.z<<"|"<<endl;

	return tmpchar.str();
}
