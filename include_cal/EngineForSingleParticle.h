//
// Created by ylang on 17-12-19.
//

#ifndef LEPLASMA_VISUAL_ENGINEFORSINGLEPARTICLE_H
#define LEPLASMA_VISUAL_ENGINEFORSINGLEPARTICLE_H

#include "Engine.h"
#include "Vector3D.h"
#include "Macros.h"
#include "Range.h"

class Grid;
class Particle;
class Range;
class Tensor3D;


class EngineForSingleParticle: public Engine{
public:
    EngineForSingleParticle(Grid*, double deltaT, Particle *uniqueParticle, map<string, double> &units);
    virtual ~EngineForSingleParticle();
    void buildCache(const Range&);

    void updateP(const Range&);
    void buildCoefficientTensor(Tensor3D& coefficientTensor,const Range &range, double deltaT, Particle *curParticle);
    void buildRHSVector(Vector3D &RHS,const Range &range, double deltaT, Particle *curParticle);
    void rootOfLinearEquationSet(Vector3D& root, const Tensor3D &coefficient, const Vector3D &RHS);

    void updateX(const Range&);

    void update(const Range&) override;

    double lightSpeed;
    Particle *uniqueParticle;

private:
    Range PtAdjacentL;
    Range PtAdjacentR;
};


#endif //LEPLASMA_VISUAL_ENGINEFORSINGLEPARTICLE_H
