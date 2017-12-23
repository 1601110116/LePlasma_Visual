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
    deltaT = 2*M_PI/(15*units["omegaCi"]);

    if(RunManager::Nodes>1){
        grid = new MPIGrid(20,20,1);
    }else{
        grid = new Grid(2,2,2,true);
    }

    particle=new Electron();

    launch(REPORT);
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
    engine=new EngineForAnalyticalField(grid,deltaT,units);

    if(RunManager::Nodes>1){
        grid->refreshParticleLocation();
    }

    if(report){
        time_t t = time(0);
        char tmp[64];
        strftime(tmp, sizeof(tmp), "%Y_%m_%d_%H_%M",localtime(&t));
        ostringstream name1;
        ostringstream name2;
        name1<<"Report"<<tmp<<"PositionX" << ".output";
        name2 << "Report" << tmp << "PositionY" << ".output";
        outFileX.open(name1.str().c_str());
        outFileY.open(name2.str().c_str());
    }
}

void CaseOfDipole::report(){
    if(REPORT){
        Particle *curParticle;
        for_each_Particle(grid,curParticle){
                            outFileX << curParticle->Position.x << '\t';
                            outFileY << curParticle->Position.y << '\t';
                        }end_for_each_Particle(curParticle)
    }
}


CaseOfDipole::~CaseOfDipole() {
    // TODO Auto-generated destructor stub
}


void CaseOfDipole::distributeParticle(){
    double keV = sqrt(2*(1.0e3*1.602e-19*1.0e7*units["gram"]*pow(units["cm"],2)*pow(units["second"],-2)));
    double x = 1.0e9*units["cm"];// * sqrt(2)/2;
    double y = 0;
    double z = 0;

    Range r=grid->World;

    Particle *p1 = particle->clone();
    p1->Position.x = x + grid->gridX()/2.0;
    p1->Position.y = y + grid->gridY()/2.0;
    p1->Position.z = z + grid->gridZ()/2.0;

    p1->Momentum.x = 10.0*keV - units["B0r03/c"]*y*pow(Square(x)+Square(y)+Square(z),-1.5);
    p1->Momentum.y = 0.0*keV + units["B0r03/c"]*x*pow(Square(x)+Square(y)+Square(z),-1.5);
    p1->Momentum.z = 10.0*keV;


    //UniformB
//    p1->Momentum.x = 10.0*keV - 1.0e-2*units["gauss"]*y/2/units["c"];
//    p1->Momentum.y = 0.0*keV + 1.0e-2*units["gauss"]*x/2/units["c"];
//    p1->Momentum.z = 10.0*keV;

    grid->directAddParticle(p1);


}


void CaseOfDipole::initP(){

}


void CaseOfDipole::initA(){


}


void CaseOfDipole::initY(){


}
