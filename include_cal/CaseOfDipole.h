//
// Created by ylang on 17-12-21.
//

#ifndef LEPLASMA_VISUAL_CASEOFDIPOLE_H
#define LEPLASMA_VISUAL_CASEOFDIPOLE_H

#include "Case.h"
#include "Particle.h"

class CaseOfDipole: public Case {
public:
    CaseOfDipole();
    virtual ~CaseOfDipole();

    void calcUnits();
    void distributeParticle();
    void initP();
    void initA();
    void initY();
    void launch(bool report);

    void report();

    Particle* particle;
    int particleCount;

    double thermalVelocity;
    double aVx,aVy,aVz;
    ofstream outFileX;
    ofstream outFileY;

private:
    double omegaCi;
    double gauss;
    double cm;
    double gram;
    double second;
    Particle *uniqueParticle;
};


#endif //LEPLASMA_VISUAL_CASEOFDIPOLE_H
