//
// Created by ylang on 17-12-19.
//

#include <Cell.h>
#include <EngineForSingleParticle.h>
#include <Grid.h>
#include <MPIGrid.h>
#include <Particle.h>
#include <Range.h>
#include <Tensor3D.h>
#include <Vertex.h>
#include <cmath>
#include <iostream>
#include "LePlasma.h"


#include <indexCache.h>




EngineForSingleParticle::EngineForSingleParticle(Grid *_grid, double dt): Engine(_grid, dt) {
#if MPI_PARALLEL && USE_CACHE
    PtAdjacentL=Range(P_ADJ_BEGIN_L,SWAP_LEN,0,grid->_height,0,grid->_length);
	PtAdjacentR=Range(grid->_width-SWAP_LEN,grid->_width-P_ADJ_BEGIN_R,0,grid->_height,0,grid->_length);
#endif

}

EngineForSingleParticle::~EngineForSingleParticle() {

}

void EngineForSingleParticle::update(const Range &range) {
#if USE_CACHE

    buildCache(range);

#if MPI_PARALLEL

#pragma omp sections
	{
#pragma omp section
		{
			buildCache(PtAdjacentL);
		}
#pragma omp section
		{
			buildCache(PtAdjacentR);
		}
	}

#endif


#pragma omp barrier //wait until all processes are finished

#endif
    updateP(range);
#pragma omp barrier

    //in MPI running, all Nodes need to sync P in shared areas here.
#if MPI_PARALLEL
    #pragma omp single
	{
		grid->syncSharedParticles();
	}
#endif
#pragma omp barrier
    updateX(range);
#pragma omp barrier

}

#define fastDecline1InLoop(i,loop) (i-1>=0?i-1:loop-1)
#define fastPlus1InLoop(i,loop) (i+1<loop?i+1:0)
void EngineForSingleParticle::rootOfLinearEquationSet(Vector3D& root, const Tensor3D &coefficient, const Vector3D &RHS){
    Tensor3D adjacentMatrix;
    double determinant = 0;

    adjacentMatrix.x.x = coefficient.y.y * coefficient.z.z - coefficient.y.z * coefficient.z.y;
    adjacentMatrix.x.y = -coefficient.x.y * coefficient.z.z + coefficient.x.z * coefficient.z.y;
    adjacentMatrix.x.z = coefficient.x.y * coefficient.y.z - coefficient.x.z * coefficient.y.y;
    adjacentMatrix.y.x = -coefficient.y.x * coefficient.z.z + coefficient.y.z * coefficient.z.x;
    adjacentMatrix.y.y = coefficient.x.x * coefficient.z.z - coefficient.x.z * coefficient.z.x;
    adjacentMatrix.y.z = -coefficient.x.x * coefficient.y.z + coefficient.x.z * coefficient.y.x;
    adjacentMatrix.z.x = coefficient.y.x * coefficient.z.y - coefficient.y.y * coefficient.z.x;
    adjacentMatrix.z.y = -coefficient.x.x * coefficient.z.y + coefficient.x.y * coefficient.z.x;
    adjacentMatrix.z.z = coefficient.x.x * coefficient.y.y - coefficient.x.y * coefficient.y.x;

    determinant += coefficient.x.x * adjacentMatrix.x.x;
    determinant += coefficient.y.x * adjacentMatrix.x.y;
    determinant += coefficient.z.x * adjacentMatrix.x.z;

    if (fabs(determinant) < 1e-16){
        std::cerr << "Coefficient Matrix is sigular, unable to find solution." << std::endl;
        std::cerr<<"Determinant: "<<determinant<<std::endl;
        //exit(0);
    }

#if OPTIMIZE
    adjacentMatrix/=determinant;

	root.x = adjacentMatrix.x * RHS;
	root.y = adjacentMatrix.y * RHS;
	root.z = adjacentMatrix.z * RHS;
#else
    Tensor3D inverseMatrix;

    inverseMatrix = adjacentMatrix/determinant;

    root.x = inverseMatrix.x * RHS;
    root.y = inverseMatrix.y * RHS;
    root.z = inverseMatrix.z * RHS;
#endif
}

void EngineForSingleParticle::buildCache(const Range& r){

    Vertex *v;
    Vector3D VertexRealPosition;
    Vector3D distance;

    Particle* p;

    for_each_Particle_within(grid,p,r){

                        for_each_Vertex_around(grid, v,p, VertexRealPosition){

                                        //d=p-v
                                        distance=p->Position;
                                        distance-=VertexRealPosition;
                                        p->W_cache[_l-sx][_m-sy][_n-sz]=_W(distance);
                                        GradW(p->GW_cache[_l-sx][_m-sy][_n-sz],distance);

                                    }end_for_each_Vertex_around

                    }end_for_each_Particle(p)

}


void EngineForSingleParticle::buidCoefficientTensor(Tensor3D& coefficientTensor,const Range &range, double deltaT, Particle *curParticle){
    coefficientTensor.restore();
    Vertex* curVertex;
    Vector3D colVector;
    Vector3D rowVector;

    int i,j,k;

    Cache_for_each_Vertex_around(grid, curVertex, curParticle,i,j,k){

                    colVector=fetch_GradW_withVertex(curParticle,i,j,k);

                    rowVector = curVertex->A;
                    rowVector*=deltaT;

                    coefficientTensor.x.x += colVector.x * rowVector.x;
                    coefficientTensor.x.y += colVector.x * rowVector.y;
                    coefficientTensor.x.z += colVector.x * rowVector.z;
                    coefficientTensor.y.x += colVector.y * rowVector.x;
                    coefficientTensor.y.y += colVector.y * rowVector.y;
                    coefficientTensor.y.z += colVector.y * rowVector.z;
                    coefficientTensor.z.x += colVector.z * rowVector.x;
                    coefficientTensor.z.y += colVector.z * rowVector.y;
                    coefficientTensor.z.z += colVector.z * rowVector.z;

                }end_for_each_Vertex_around
    coefficientTensor.x.x += 1;
    coefficientTensor.y.y += 1;
    coefficientTensor.z.z += 1;
}

void EngineForSingleParticle::buildRHSVector(Vector3D &RHS,const Range &range, double deltaT, Particle *curParticle){
    RHS.restore();
    Vertex *curVertex1;
    Vertex *curVertex2;
    Vector3D dRHS;

    int i,j,k;
    int l,m,n;

    Cache_for_each_Vertex_around(grid, curVertex1, curParticle,i,j,k){

                    Cache_for_each_Vertex_around(grid, curVertex2, curParticle,l,m,n){

                                    dRHS=fetch_GradW_withVertex(curParticle,l,m,n);
                                    dRHS*=fetch_W_withVertex(curParticle,i,j,k)* (curVertex1->A * curVertex2->A);
                                    RHS +=  dRHS;

                                    //RHS +=GradW(curParticle->Position - Vertex2RealPosition) * (_W(curParticle->Position - Vertex1RealPosition) * (curVertex1->A * curVertex2->A));

                                }end_for_each_Vertex_around

                }end_for_each_Vertex_around

    RHS *= -deltaT;

    RHS +=curParticle->Momentum;
}


void EngineForSingleParticle::updateP(const Range& range){
    Particle *curParticle;
    Tensor3D coefficient;    // The 3*3 coefficient matrix
    Vector3D RHS;                // Right hand side

    for_each_Particle_within(grid, curParticle, range){
                        buidCoefficientTensor(coefficient, range, deltaT, curParticle);
                        buildRHSVector(RHS, range, deltaT, curParticle);
                        rootOfLinearEquationSet(curParticle->Momentum, coefficient, RHS);
                    }end_for_each_Particle(curParticle)

}

void EngineForSingleParticle::updateX(const Range& range){
    Particle *curParticle;

    Vertex *curVertex;
    Vector3D secondTerm;
    Vector3D dSecondTerm,dX;

    int i,j,k;

    for_each_Particle_within(grid, curParticle, range){
                        secondTerm.restore();

                        Cache_for_each_Vertex_around(grid, curVertex,curParticle, i,j,k){

                                        dSecondTerm=curVertex->A;dSecondTerm*=fetch_W_withVertex(curParticle,i,j,k);
                                        secondTerm += dSecondTerm;

                                        //secondTerm += curVertex->A * _W(curParticle->Position - VertexRealPosition);

                                    }end_for_each_Vertex_around

                        dX=curParticle->Momentum;
                        dX+=secondTerm;
                        dX*=deltaT;

                        curParticle->Position += dX;

                        //curParticle->Position += (curParticle->Momentum + secondTerm) * deltaT;

                    }end_for_each_Particle(curParticle)

}