/*
 * Vector3D.h
 *
 *  Created on: Mar 3, 2016
 *      Author: zlstudio
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include <iostream>
#include<string>
#include<sstream>
#include<math.h>

#define init_3D_Vector_Field(field,x,y,z)\
	field = new Vector3D***[x];\
	for (int i = 0; i < x; i++){\
		field[i] = new Vector3D**[y];\
		for (int j = 0; j < y; j++){\
			field[i][j] = new Vector3D*[z];\
			for (int k = 0; k < z; k++){\
				field[i][j][k] = new Vector3D();\
			}\
		}\
	}\

#define init_3D_Scalar_Field(field,x,y,z)\
	field = new double**[x];\
	for (int i = 0; i < x; i++){\
		field[i] = new double*[y];\
		for (int j = 0; j < y; j++){\
			field[i][j] = new double[z];\
			for (int k = 0; k < z; k++){\
				field[i][j][k] = 0;\
			}\
		}\
	}\

using namespace std;

class IV_Vector3D;

class Vector3D {
	public:
		Vector3D();
		Vector3D(double x,double y,double z);
		virtual ~Vector3D();
		double x,y,z;

		///////ANTI-MEMORY LEAK & FAST IV CALCULATE
		IV_Vector3D operator +(const Vector3D &a)const;
		IV_Vector3D operator -(const Vector3D &a)const;
		IV_Vector3D operator *(const double &a)const;
		IV_Vector3D operator /(const double &a)const;


		//SAFE
		Vector3D &operator =(const Vector3D &a);

		Vector3D &operator +=(const Vector3D &a);

		Vector3D &operator -=(const Vector3D &a);

		double operator *(const Vector3D &a)const;

		double operator >>(const Vector3D &a)const;

		//Very safe
		Vector3D &operator *=(const double &a);
		Vector3D &operator /=(const double &a);

		void restore();

		string toString();
};

class IV_Vector3D: public Vector3D{
	public:
		IV_Vector3D();
		IV_Vector3D(double x,double y,double z);
		virtual ~IV_Vector3D();

		IV_Vector3D &operator +(const Vector3D &a);
		IV_Vector3D &operator -(const Vector3D &a);

		IV_Vector3D &operator *(const double &a);
		IV_Vector3D &operator /(const double &a);

};

#endif /* VECTOR3D_H_ */
