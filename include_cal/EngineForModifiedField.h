//
// Created by ylang on 17-12-28.
//

#ifndef LEPLASMA_VISUAL_ENGINEFORMODIFIEDFIELD_H
#define LEPLASMA_VISUAL_ENGINEFORMODIFIEDFIELD_H

#include "Engine.h"
#include "Tensor3D.h"
#include "Particle.h"
#include <map>

class EngineForModifiedField: public Engine {
public:
    EngineForModifiedField(Grid *grid, double deltaT, map<string, double> &units);
    virtual ~EngineForModifiedField();

//    double x,y,z;  // used to calc local field of the particle
    double r2, r_3, ax2, ay2, r_5, r0_1;
    static Tensor3D gradAdt;
    void calcA(Particle *p);
    void calcGradAdt(Particle *p);
    void updateP(const Range&);
    void updateX(const Range&);
    void buildCoefficientTensor(Tensor3D& coefficientTensor);
    void rootOfLinearEquationSet(Vector3D &root, const Tensor3D &coefficient, const Vector3D &RHS);

    void update(const Range&) override ;

};


#endif //LEPLASMA_VISUAL_ENGINEFORMODIFIEDFIELD_H
