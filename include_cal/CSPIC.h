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
#include "Macros.h"
#include "Range.h"

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
    double lightSpeed;

    CSPIC(Grid*,double dt, double lightSpeed);
    virtual ~CSPIC();

private:
    Range PtAdjacentL;
    Range PtAdjacentR;
};




#endif /* INCLUDE_CSPIC_H_ */
