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

inline double _W(const Vector3D &r){
    return _W1(r.x)*_W1(r.y)*_W1(r.z);
}


void CaseOfDipole::calcUnits() {
    // cell width in cm
    double dx = 4.0e7;
    // super particle weight in g
    double mp = 1.6726e-24;
    // super particle charge in g^1/2cm^3/2s^-1
    double qp = 4.8032e-10;
    cm = 1.0/dx;
    gram = 1.0/mp;
    // second in mp^1/2dx^3/2qp^-1
    second = qp*pow(gram,0.5)*pow(cm,1.5);
    lightSpeed = 2.9979e10*cm/second;
    cout << "lightSpeed = " << lightSpeed << endl;

    gauss = pow(gram,0.5) * pow(cm,-0.5) * pow(second,-1);
    omegaCi = 9.58e3*1.0e-2/second;

}



CaseOfDipole::CaseOfDipole() {
    calcUnits();

/*	You should disable "updateA(range)" and "updateY(range)" first */
    deltaT = 1/(8*omegaCi);

    if(RunManager::Nodes>1){
        grid = new MPIGrid(20,20,1);
    }else{
        grid = new Grid(64,64,64,true);
    }

    particle=new Electron();
    particleCount=1;

    thermalVelocity=0;//0.1*LIGHT_SPEED;

    double keV = sqrt(2*(1.0e3*1.602e-19*1.0e7*gram*pow(cm,2)*pow(second,-2)));
    aVx=0;
    aVy=10*keV;
    aVz = 10*keV;

    //select Engine

    engine=new EngineForSingleParticle(grid,deltaT,lightSpeed);

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


inline double Random(double min, double max){
    return (rand()*(max-min))/RAND_MAX+min;
}

static double gaussRand(double expectation, double sigma){
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if (phase == 0) {
        do{
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        }  while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else
        X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    X = X * sigma + expectation;
    return X;
}

void CaseOfDipole::distributeParticle(){

    double sigma = thermalVelocity / sqrt(2);

    Range r=grid->World;

    for (int i = 0; i < particleCount; i++){
        Particle *newParticle = particle->clone();

        newParticle->Position.x=1e9*cm + grid->gridX()/2.0;
        newParticle->Position.y=0.0 + grid->gridY()/2.0;
        newParticle->Position.z=0.0 + grid->gridZ()/2.0;

        newParticle->Momentum.x=gaussRand(aVx, sigma);
        newParticle->Momentum.y=gaussRand(aVy, sigma);
        newParticle->Momentum.z=gaussRand(aVz, sigma);

        grid->directAddParticle(newParticle);
    }

}


void CaseOfDipole::initP(){
    Particle *curParticle;
    Vertex *curVertex;

    //Vector3D v1, v2,
    Vector3D VertexRealPosition;
    for_each_Particle_within(grid, curParticle,grid->workSpace){
                        //currently Momentum is Velocity
                        for_each_Vertex_around(grid, curVertex, curParticle,VertexRealPosition){
                                        curParticle->Momentum += curVertex->A * _W(curParticle->Position - VertexRealPosition);
                                    }end_for_each_Vertex_around

                    }end_for_each_Particle(curParticle)
}


void CaseOfDipole::initA(){
    //	Particle *curParticle;
    //	double rc2;  //the distance between the source point and the field point multimplied by c square
    //
    //	Vertex *curVertex;
    //	Vector3D VertexRealPosition;
    //	Vector3D r;
    //	for_each_Particle_within(grid, curParticle,grid->workSpace){
    //		//currently Momentum is Velocity
    //		for_each_Vertex_around(grid, curVertex, curParticle, VertexRealPosition){
    //			r = VertexRealPosition - curParticle->Position;
    //			rc2 = sqrt(Square(r.x) + Square(r.y) + Square(r.z)) * Square(LIGHT_SPEED);
    //			curVertex->A += (curParticle->Momentum / rc2);
    //		}end_for_each_Vertex_around
    //
    //	}end_for_each_Particle(curParticle)
    Vertex *curVertex;
    double B0r03 = 1.0e-2*gauss * pow(1.0e9*cm,3);
    double Ax, Ay, Az;
    double x,y,z;
    double r,theta,phi;
//    double B0r03 = 1.0e-2*gauss*Cube(1e9);
    for_each_Vertex_within(grid,curVertex,grid->workSpace){
                    x = curVertex->x()-grid->gridX()/2.0;
                    y = curVertex->y()-grid->gridY()/2.0;
                    z = curVertex->z()-grid->gridZ()/2.0;
                    r = sqrt(Square(x)+Square(y)+Square(z));
                    theta = acos(z/r);
                    phi = atan(y/x);
                    Ax = -sin(theta)*sin(phi);
                    Ay = sin(theta)*cos(phi);
                    Az = 0.0;
                    curVertex->A=Vector3D(Ax, Ay, Az)*(B0r03/(Square(r)*lightSpeed));
                }end_for_each_Vertex_within
}


void CaseOfDipole::initY(){
    //	Particle *curParticle;
    //	Vector3D r;  //pointing from source point to field point
    //	double lengthR;  //the length of r
    //
    //
    //	Vertex *curVertex;
    //	Vector3D VertexRealPosition;
    //	for_each_Particle_within(grid, curParticle, grid->workSpace){
    //		for_each_Vertex_around(grid, curVertex, curParticle, VertexRealPosition){
    //			r = VertexRealPosition - curParticle->Position;
    //			lengthR = sqrt(Square(r.x) + Square(r.y) + Square(r.z));
    //			curVertex->Y = r / (-Cube(lengthR));
    //		}end_for_each_Vertex_around
    //	}end_for_each_Particle(curParticle)


}
