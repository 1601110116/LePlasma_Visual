//
// Created by ylang on 17-12-22.
//

#ifndef LEPLASMA_VISUAL_ENGINEFORANALYTICALFIELD_H
#define LEPLASMA_VISUAL_ENGINEFORANALYTICALFIELD_H

#include "Engine.h"
#include "Tensor3D.h"
#include "Particle.h"
#include <map>

class EngineForAnalyticalField: public Engine {
public:
    EngineForAnalyticalField(Grid *grid, double deltaT, map<string, double> &units);
    virtual ~EngineForAnalyticalField();

    double x,y,z;  // used to calc local field of the particle
    double tmp1, tmp2, tmp3;
    Tensor3D gradAdt;
    void calcA(Particle *p);
    void calcGradAdt();
    void updateP(const Range&);
    void updateX(const Range&);
    void buildCoefficientTensor(Tensor3D& coefficientTensor);
    void rootOfLinearEquationSet(Vector3D &root, const Tensor3D &coefficient, const Vector3D &RHS);

    void update(const Range&) override ;

};

#endif //LEPLASMA_VISUAL_ENGINEFORANALYTICALFIELD_H
