//
// Created by ylang on 17-12-21.
//


#include <Cell.h>
#include <CaseOfDipole.h>
#include <Electron.h>
#include <Grid.h>
#include <LePlasma.h>
#include <MPIGrid.h>
#include <Macros.h>
#include <RunManager.h>
#include <stdlib.h>
#include <Vector3D.h>
#include <Vertex.h>
#include <cmath>
#include <iostream>
#include "EngineForSingleParticle.h"
#include "EngineForAnalyticalField.h"
#include "EngineForModifiedField.h"


void CaseOfDipole::calcUnits() {
    // cell width in cm
    double dx = 2.0e9;
    // super particle weight in g
    double mp = 1.6726e-24;
    // super particle charge in g^1/2cm^3/2s^-1
    double qp = 4.8032e-10;

    units["cm"] = 1.0/dx;
    units["gram"] = 1.0/mp;
    units["second"] = qp*pow(units["gram"],0.5)*pow(units["cm"],1.5);
    units["c"] = 2.9979e10*units["cm"]/units["second"];
    units["gauss"] = pow(units["gram"],0.5) * pow(units["cm"],-0.5) * pow(units["second"],-1);
    units["omegaCi"] = 9.58e3*1.0e-2/units["second"];

    units["B0r03/c"] = 1.0e-2*units["gauss"] * pow(1.0e9*units["cm"],3) / units["c"];
    cout << "lightSpeed = " << units["c"] << endl;
    cout << "B0r0 = " << 1e-2*units["gauss"]*1e9*units["cm"] << endl;
}

CaseOfDipole::CaseOfDipole() {
    calcUnits();
    deltaT = 2*M_PI/(40*units["omegaCi"]);

    if(RunManager::Nodes>1){
        grid = new MPIGrid(20,20,1);
    }else{
        grid = new Grid(2,2,2,true);
    }

    particle=new Electron();

//    launch(REPORT);
    launch(true);
}

void CaseOfDipole::launch(bool report){

    srand(0);

    if(RunManager::MPI_ID==0){
        distributeParticle();
    }

    grid->refreshParticleLocation();

    cout<<"Node "<<RunManager::MPI_ID<<" Distribute Complete. Particles: "<<grid->particles()<<endl;

    initA();
    initP();
    initY();

    //select Engine
//    engine=new EngineForModifiedField(grid,deltaT,units);
    engine=new EngineForModifiedField(grid,deltaT,units);

    if(RunManager::Nodes>1){
        grid->refreshParticleLocation();
    }

    if(report){
        time_t t = time(0);
        char tmp[64];
        strftime(tmp, sizeof(tmp), "%Y_%m_%d_%H_%M",localtime(&t));
        ostringstream name1;
        ostringstream name2;
        name1<<"Report"<<tmp<<"Position" << ".output";
        name2 << "Report" << tmp << "Velocity" << ".output";
        outFileX.open(name1.str().c_str());
        outFileV.open(name2.str().c_str());
    }
}

void CaseOfDipole::report(){
//    if(REPORT){
        Particle *p;
        for_each_Particle(grid,p){
                            outFileX << p->X.x/units["cm"] << '\t' << p->X.y/units["cm"] << '\t' << p->X.z/units["cm"] << endl;
                            outFileV << (p->Momentum.x-p->A.x)*units["second"]/units["cm"] << '\t';
                            outFileV << (p->Momentum.y-p->A.y)*units["second"]/units["cm"] << '\t';
                            outFileV << (p->Momentum.z-p->A.z)*units["second"]/units["cm"] << endl;

//                            outFileX << curParticle->Position.x << '\t';
//                            outFileY << curParticle->Position.y << '\t';
                        }end_for_each_Particle(p)
//    }
}


CaseOfDipole::~CaseOfDipole() {
    // TODO Auto-generated destructor stub
}


void CaseOfDipole::distributeParticle(){
    double v0 = sqrt(2*(1.0e4*1.602e-19*1.0e7*units["gram"]*pow(units["cm"],2)*pow(units["second"],-2)));
    Range r=grid->World;

    Particle *p1 = particle->clone();
    p1->X.x = 1.0e9*units["cm"];// * sqrt(2)/2;
    p1->X.y = 0;
    p1->X.z = 0;
    p1->Position.x = p1->X.x + grid->gridX()/2.0;
    p1->Position.y = p1->X.y + grid->gridY()/2.0;
    p1->Position.z = p1->X.z + grid->gridZ()/2.0;

    p1->A.x = -1*units["B0r03/c"]*p1->X.y*pow(Square(p1->X.x)+Square(p1->X.y)+Square(p1->X.z),-1.5);
    p1->A.y = units["B0r03/c"]*p1->X.x*pow(Square(p1->X.x)+Square(p1->X.y)+Square(p1->X.z),-1.5);
    p1->A.z = 0;
    p1->Momentum.x = v0 + p1->A.x;
    p1->Momentum.y = 0.0 + p1->A.y;
    p1->Momentum.z = v0 + p1->A.z;


    //UniformB
//    p1->Momentum.x = 10.0*keV - 1.0e-2*units["gauss"]*y/2/units["c"];
//    p1->Momentum.y = 0.0*keV + 1.0e-2*units["gauss"]*x/2/units["c"];
//    p1->Momentum.z = 10.0*keV;

    grid->directAddParticle(p1);


}


void CaseOfDipole::initP(){

}


void CaseOfDipole::initA() {

}


void CaseOfDipole::initY(){


}
