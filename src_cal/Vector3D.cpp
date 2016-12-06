/*
 * Vector3D.cpp
 *
 *  Created on: Mar 3, 2016
 *      Author: zlstudio
 */

#include "Vector3D.h"

Vector3D::Vector3D():x(0),y(0),z(0) {

}
Vector3D::Vector3D(double x,double y,double z):x(x),y(y),z(z) {

}

Vector3D::~Vector3D() {
}

/////Intermediate Variable

IV_Vector3D::IV_Vector3D(){
}

IV_Vector3D::IV_Vector3D(double x,double y,double z):Vector3D(x,y,z){
}

IV_Vector3D::~IV_Vector3D() {
}

///////ANTI-MEMORY LEAK & FAST IV CALCULATE

IV_Vector3D Vector3D::operator +(const Vector3D &a)const{
	return IV_Vector3D(a.x+x,a.y+y,a.z+z);
}

IV_Vector3D Vector3D::operator -(const Vector3D &a)const{
	return IV_Vector3D(x-a.x,y-a.y,z-a.z);
}

IV_Vector3D Vector3D::operator *(const double &a)const{
	return IV_Vector3D(x*a,y*a,z*a);
};

IV_Vector3D Vector3D::operator /(const double &a)const{
	return IV_Vector3D(x/a,y/a,z/a);
}

///////////////////////////////////////////////IV (return this)

IV_Vector3D &IV_Vector3D::operator +(const Vector3D &a){
	x+=a.x;
	y+=a.y;
	z+=a.z;
	return *this;
}

IV_Vector3D &IV_Vector3D::operator -(const Vector3D &a){
	x-=a.x;
	y-=a.y;
	z-=a.z;
	return *this;
}

IV_Vector3D &IV_Vector3D::operator *(const double &a){
	x*=a;
	y*=a;
	z*=a;
	return *this;
};

IV_Vector3D &IV_Vector3D::operator /(const double &a){
	x/=a;
	y/=a;
	z/=a;
	return *this;
}


/////Safe (return this)
Vector3D &Vector3D::operator +=(const Vector3D &a){
	x+=a.x;
	y+=a.y;
	z+=a.z;
	return *this;
}

Vector3D &Vector3D::operator -=(const Vector3D &a){
	x-=a.x;
	y-=a.y;
	z-=a.z;
	return *this;
}

Vector3D &Vector3D::operator *=(const double &a){

	x*=a;
	y*=a;
	z*=a;
	return *this;
}

Vector3D &Vector3D::operator /=(const double &a){
	x/=a;
	y/=a;
	z/=a;
	return *this;
}

Vector3D &Vector3D::operator =(const Vector3D &a){
	x=a.x;
	y=a.y;
	z=a.z;
	return *this;
}


//Safe (return double)
double Vector3D::operator *(const Vector3D &a)const{
	return a.x*x+a.y*y+a.z*z;
}


double Vector3D::operator >>(const Vector3D &a)const {
	return sqrt((x-a.x)*(x-a.x)+(y-a.y)*(y-a.y)+(z-a.z)*(z-a.z));
}

void Vector3D::restore(){
	x=y=z=0;
}
string Vector3D::toString(){

	ostringstream tmpchar;

	tmpchar<<"("<<x<<","<<y<<","<<z<<")";

	return tmpchar.str();
}
