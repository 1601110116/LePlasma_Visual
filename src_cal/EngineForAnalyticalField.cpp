//
// Created by ylang on 17-12-22.
//

//　This Engine requires OMP and MPI parallels to be turned off

#include "EngineForAnalyticalField.h"
#include "Grid.h"
#include "CaseOfDipole.h"

EngineForAnalyticalField::EngineForAnalyticalField(Grid *grid, double deltaT, map<string, double> &units):\
 Engine(grid, deltaT, units){

}

EngineForAnalyticalField::~EngineForAnalyticalField() {

}

void EngineForAnalyticalField::update(const Range &range){
    updateP(range);
    updateX(range);
}

void EngineForAnalyticalField::updateP(const Range &range) {
    Tensor3D coefficient;
    Vector3D RHS;
    Particle *p;
    for_each_Particle_within(grid, p, range){
                        buildCoefficientTensor(coefficient, p);
                        buildRHSVector(RHS, p);
                        rootOfLinearEquationSet(p->Momentum, coefficient, RHS);
                    }end_for_each_Particle(p)
}

void EngineForAnalyticalField::updateX(const Range &range) {

}

void EngineForAnalyticalField::buildCoefficientTensor(Tensor3D &coefficientTensor, Particle *curParticle) {
    double x = curParticle->Position.x - grid->gridX()/2.0;
    double y = curParticle->Position.y - grid->gridY()/2.0;
    double z = curParticle->Position.z - grid->gridZ()/2.0;
    coefficientTensor.restore();
    coefficientTensor.x.x =1;
}

void EngineForAnalyticalField::buildRHSVector(Vector3D &RHS, Particle *curParticle) {

}

void EngineForAnalyticalField::rootOfLinearEquationSet(Vector3D &root, const Tensor3D &coefficient,
                                                       const Vector3D &RHS) {

}