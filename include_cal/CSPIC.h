/*
 * CSPIC.h
 *
 *  Created on: 2016-4-18
 *      Author: lyt
 */

#ifndef INCLUDE_CSPIC_H_
#define INCLUDE_CSPIC_H_


#include "Engine.h"
#include "Vector3D.h"
#include "PhysicalParameters.h"

class Grid;
class Particle;
class Range;
class Tensor3D;



class CSPIC:public Engine{
	public:

		Vector3D ****curldA;
		Vector3D ****curldTCurldA;

		void buildCache(const Range&);


		void updateP(const Range&);
		void buidCoefficientTensor(Tensor3D& coefficientTensor,const Range &range, double deltaT, Particle *curParticle);
		void buildRHSVector(Vector3D &RHS,const Range &range, double deltaT, Particle *curParticle);
		void rootOfLinearEquationSet(Vector3D& root, const Tensor3D &coefficient, const Vector3D &RHS);


		void updateY(const Range&);
		void curld(Vector3D& result, Vector3D ****vectorField, int i, int j, int k);
		void curldT(Vector3D& result,Vector3D ****vectorField, int i, int j, int k);
		void setCurldA(const Range&);
		void setCurldTCurldA(const Range&);

		void updateX(const Range&);
		void updateA(const Range&);
		void updatePFalse(const Range&);

		void update(const Range&);

		CSPIC(Grid*,double);
		virtual ~CSPIC();

};


inline double _W1(double x){
	if (x > 2)
		return 0.0;
	else if (x > 1)
		return x * (Cube(x) * (x * (x * (x * (15.0/1024 * x - 15.0/128) + 49.0/128) - 21.0/32) + 35.0/64) - 1.0) + 1.0;
	else if (x > 0)
		return Square(x) * (Square(x) * (x * (x * (x * (-15.0/1024 * x - 15.0/128) + 7.0/16) - 21.0/32) + 175.0/256) - 105.0/128) + 337.0/512;
	else if (x > -1)
		return x * x * (x * x * (x * (x * (x * (-15.0/1024 * x + 15.0/128) + 7.0/16) + 21.0/32) + 175.0/256) - 105.0/128) + 337.0/512;
	else if (x > -2)
		return x * (Cube(x) * (x * (x * (x * (15.0/1024 * x + 15.0/128) + 49.0/128) + 21.0/32) + 35.0/64) + 1.0) + 1.0;
	else
		return 0.0;
}

inline double dW1(double x){
	if (x > 2)
		return 0.0;
	else if (x > 1)
		return Cube(x) * (x * (x * (x * (15.0/128 * x - 105.0/128) + 147.0/64) - 105.0/32) + 35.0/16) - 1.0;
	else if (x > 0)
		return x * (x * x * (x * (x * (x * (-15.0/128 * x - 105.0/128) + 21.0/8) - 105.0/32) + 175.0/64) - 105.0/64);
	else if (x > -1)
		return x * (x * x * (x * (x * (x * (-15.0/128 * x + 105.0/128) + 21.0/8) + 105.0/32) + 175.0/64) - 105.0/64);
	else if (x > -2)
		return Cube(x) * (x * (x * (x * (15.0/128 * x + 105.0/128) + 147.0/64) + 105.0/32) + 35.0/16) + 1.0;
	else
		return 0.0;
}



inline double _W(const Vector3D &r){
	return _W1(r.x)*_W1(r.y)*_W1(r.z);
}


inline void GradW(Vector3D& gradW, const Vector3D &r){
	gradW.x=dW1(r.x)*_W1(r.y)*_W1(r.z);
	gradW.y=_W1(r.x)*dW1(r.y)*_W1(r.z);
	gradW.z=_W1(r.x)*_W1(r.y)*dW1(r.z);
}

#endif /* INCLUDE_CSPIC_H_ */
