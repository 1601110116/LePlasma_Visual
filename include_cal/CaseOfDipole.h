//
// Created by ylang on 17-12-21.
//

#ifndef LEPLASMA_VISUAL_CASEOFDIPOLE_H
#define LEPLASMA_VISUAL_CASEOFDIPOLE_H

#include "Case.h"
#include "Particle.h"
#include <map>
#include <string>

class CaseOfDipole: public Case {
public:
    CaseOfDipole();
    virtual ~CaseOfDipole();

    void calcUnits() override ;
    void distributeParticle() override;
    void initP() override ;
    void initA() override ;
    void initY() override ;
    void launch(bool report) override;

    void report() override;

    Particle* particle;
    int particleCount;

    double thermalVelocity;
    double aVx,aVy,aVz;
    ofstream outFileX;
    ofstream outFileY;


    map<string, double> units;
    double omegaCi;
    double gauss;
    double cm;
    double gram;
    double second;
    Particle *uniqueParticle;
};


#endif //LEPLASMA_VISUAL_CASEOFDIPOLE_H
