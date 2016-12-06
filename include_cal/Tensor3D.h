/*
 * Tensor3D.h
 *
 *  Created on: Mar 22, 2016
 *      Author: zlstudio
 */

#ifndef TENSOR3D_H_
#define TENSOR3D_H_

#include "Vector3D.h"

#define init_3D_Tensor_Field(field,x,y,z)\
	field = new Tensor3D***[x];\
	for (int i = 0; i < x; i++){\
		field[i] = new Tensor3D**[y];\
		for (int j = 0; j < y; j++){\
			field[i][j] = new Tensor3D*[z];\
			for (int k = 0; k < z; k++){\
				field[i][j][k] = new Tensor3D();\
			}\
		}\
	}\

class Tensor3D {
	public:
		Tensor3D();
		Tensor3D(Vector3D x, Vector3D y, Vector3D z);
		virtual ~Tensor3D();

		Vector3D x,y,z;

		Tensor3D &operator =(const Tensor3D &a);
		Tensor3D &operator *=(double a);
		Tensor3D &operator /=(double a);

		Tensor3D operator /(double a)const;

		void restore();

		string toString();

};

#endif /* TENSOR3D_H_ */
