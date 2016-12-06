/*
 * CSPIC.cpp
 *
 *  Created on: 2016-4-18
 *      Author: lyt
 */

#include <Cell.h>
#include <CSPIC.h>
#include <Grid.h>
#include <MPIGrid.h>
#include <Particle.h>
#include <Range.h>
#include <Tensor3D.h>
#include <Vertex.h>
#include <cmath>
#include <iostream>

#if USE_CACHE
#include <indexCache.h>
#endif

Range PtAdjacentL,PtAdjacentR;

CSPIC::CSPIC(Grid*_grid,double dt):Engine(_grid,dt){

	init_3D_Vector_Field(curldA,grid->_width,grid->_height,grid->_length);
	init_3D_Vector_Field(curldTCurldA, grid->_width,grid->_height,grid->_length);
	std::cout << "CSPIC initialized" << std::endl;

#if MPI_PARALLEL && USE_CACHE
	PtAdjacentL=Range(P_ADJ_BEGIN_L,SWAP_LEN,0,grid->_height,0,grid->_length);
	PtAdjacentR=Range(grid->_width-SWAP_LEN,grid->_width-P_ADJ_BEGIN_R,0,grid->_height,0,grid->_length);
#endif

}

CSPIC::~CSPIC(){

}

void CSPIC::update(const Range& range){

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


	//start update
	updateP(range);
//	updatePFalse(range);
#pragma omp barrier

	//in MPI running, all Nodes need to sync P in shared areas here.
#if MPI_PARALLEL
#pragma omp single
	{
		grid->syncSharedParticles();
	}
#endif


#pragma omp barrier
	setCurldA(range);
#pragma omp barrier
	setCurldTCurldA(range);

	updateY(range);


#pragma omp barrier
	updateX(range);


#pragma omp barrier
	updateA(range);
	//implicit omp barrier
}

void CSPIC::setCurldA(const Range& range){
	int i, j, k;
	in_Range(i,j,k,range){
		CSPIC::curld(*curldA[i][j][k], grid->A_indexer, i, j, k);
	}end_in_Range
}

void CSPIC::setCurldTCurldA(const Range& range){
	int i, j, k;
	in_Range(i,j,k,range){
		CSPIC::curldT(*curldTCurldA[i][j][k],curldA, i, j, k);
	}end_in_Range
}

#define fastDecline1InLoop(i,loop) (i-1>=0?i-1:loop-1)
#define fastPlus1InLoop(i,loop) (i+1<loop?i+1:0)

void CSPIC::curld(Vector3D& result, Vector3D ****vectorField, int i, int j, int k){

#if OPTIMIZE_1D_MPI
	result.x = (vectorField[i][j][k]->z - vectorField[i][0][k]->z)-(vectorField[i][j][k]->y - vectorField[i][j][0]->y);
	result.y = (vectorField[i][j][k]->x - vectorField[i][j][0]->x)-(vectorField[i][j][k]->z - vectorField[i-1][j][k]->z);
	result.z = (vectorField[i][j][k]->y - vectorField[i-1][j][k]->y)-(vectorField[i][j][k]->x - vectorField[i][0][k]->x);
#elif OPTIMIZE_1D
	result.x = (vectorField[i][j][k]->z - vectorField[i][0][k]->z)-(vectorField[i][j][k]->y - vectorField[i][j][0]->y);
	result.y = (vectorField[i][j][k]->x - vectorField[i][j][0]->x)-(vectorField[i][j][k]->z - vectorField[fastDecline1InLoop(i,grid->_width)][j][k]->z);
	result.z = (vectorField[i][j][k]->y - vectorField[fastDecline1InLoop(i,grid->_width)][j][k]->y)-(vectorField[i][j][k]->x - vectorField[i][0][k]->x);
#elif OPTIMIZE_3D
	result.x = (vectorField[i][j][k]->z - vectorField[i][fastDecline1InLoop(j,grid->_height)][k]->z)-(vectorField[i][j][k]->y - vectorField[i][j][fastDecline1InLoop(k,grid->_length)]->y);
	result.y = (vectorField[i][j][k]->x - vectorField[i][j][fastDecline1InLoop(k,grid->_length)]->x)-(vectorField[i][j][k]->z - vectorField[fastDecline1InLoop(i,grid->_width)][j][k]->z);
	result.z = (vectorField[i][j][k]->y - vectorField[fastDecline1InLoop(i,grid->_width)][j][k]->y)-(vectorField[i][j][k]->x - vectorField[i][fastDecline1InLoop(j,grid->_height)][k]->x);
#else
	result.x = (vectorField[i][j][k]->z - vectorField[i][(j-1+grid->_height)%grid->_height][k]->z)-(vectorField[i][j][k]->y - vectorField[i][j][(k-1+grid->_length)%grid->_length]->y);
	result.y = (vectorField[i][j][k]->x - vectorField[i][j][(k-1+grid->_length)%grid->_length]->x)-(vectorField[i][j][k]->z - vectorField[(i-1+grid->_width)%grid->_width][j][k]->z);
	result.z = (vectorField[i][j][k]->y - vectorField[(i-1+grid->_width)%grid->_width][j][k]->y)-(vectorField[i][j][k]->x - vectorField[i][(j-1+grid->_height)%grid->_height][k]->x);
#endif
}

void CSPIC::curldT(Vector3D& result,Vector3D ****vectorField, int i, int j, int k){
#if OPTIMIZE_1D_MPI
	result.x = (vectorField[i][0][k]->z - vectorField[i][j][k]->z)-(vectorField[i][j][0]->y - vectorField[i][j][k]->y);
	result.y = (vectorField[i][j][0]->x - vectorField[i][j][k]->x)-(vectorField[i+1][j][k]->z - vectorField[i][j][k]->z);
	result.z = (vectorField[i+1][j][k]->y - vectorField[i][j][k]->y)-(vectorField[i][0][k]->x - vectorField[i][j][k]->x);
#elif OPTIMIZE_1D
	result.x = (vectorField[i][0][k]->z - vectorField[i][j][k]->z)-(vectorField[i][j][0]->y - vectorField[i][j][k]->y);
	result.y = (vectorField[i][j][0]->x - vectorField[i][j][k]->x)-(vectorField[fastPlus1InLoop(i,grid->_width)][j][k]->z - vectorField[i][j][k]->z);
	result.z = (vectorField[fastPlus1InLoop(i,grid->_width)][j][k]->y - vectorField[i][j][k]->y)-(vectorField[i][0][k]->x - vectorField[i][j][k]->x);
#elif OPTIMIZE_3D
	result.x = (vectorField[i][fastPlus1InLoop(j,grid->_height)][k]->z - vectorField[i][j][k]->z)-(vectorField[i][j][fastPlus1InLoop(k,grid->_length)]->y - vectorField[i][j][k]->y);
	result.y = (vectorField[i][j][fastPlus1InLoop(k,grid->_length)]->x - vectorField[i][j][k]->x)-(vectorField[fastPlus1InLoop(i,grid->_width)][j][k]->z - vectorField[i][j][k]->z);
	result.z = (vectorField[fastPlus1InLoop(i,grid->_width)][j][k]->y - vectorField[i][j][k]->y)-(vectorField[i][fastPlus1InLoop(j,grid->_height)][k]->x - vectorField[i][j][k]->x);
#else
	result.x = (vectorField[i][(j+1)%grid->_height][k]->z - vectorField[i][j][k]->z)-(vectorField[i][j][(k+1)%grid->_length]->y - vectorField[i][j][k]->y);
	result.y = (vectorField[i][j][(k+1)%grid->_length]->x - vectorField[i][j][k]->x)-(vectorField[(i+1)%grid->_width][j][k]->z - vectorField[i][j][k]->z);
	result.z = (vectorField[(i+1)%grid->_width][j][k]->y - vectorField[i][j][k]->y)-(vectorField[i][(j+1)%grid->_height][k]->x - vectorField[i][j][k]->x);
#endif
}


void CSPIC::rootOfLinearEquationSet(Vector3D& root, const Tensor3D &coefficient, const Vector3D &RHS){
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


#if USE_CACHE

void CSPIC::buildCache(const Range& r){

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


void CSPIC::buidCoefficientTensor(Tensor3D& coefficientTensor,const Range &range, double deltaT, Particle *curParticle){
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

void CSPIC::buildRHSVector(Vector3D &RHS,const Range &range, double deltaT, Particle *curParticle){
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


void CSPIC::updateP(const Range& range){
	Particle *curParticle;
	Tensor3D coefficient;    // The 3*3 coefficient matrix
	Vector3D RHS;                // Right hand side

	for_each_Particle_within(grid, curParticle, range){
		buidCoefficientTensor(coefficient, range, deltaT, curParticle);
		buildRHSVector(RHS, range, deltaT, curParticle);
		rootOfLinearEquationSet(curParticle->Momentum, coefficient, RHS);
	}end_for_each_Particle(curParticle)

}

void CSPIC::updatePFalse(const Range &range){
	Particle *curParticle;
	Vector3D firstTerm, secondTerm;
	Vertex *curVertex1, *curVertex2;
	int i,j,k;
	int l,m,n;

	for_each_Particle_within(grid, curParticle,range){
		firstTerm.restore();
		secondTerm.restore();
		Cache_for_each_Vertex_around(grid, curVertex1,curParticle,i,j,k){
			firstTerm -= fetch_GradW_withVertex(curParticle,i,j,k)*(curParticle->Momentum*curVertex1->A);
			Cache_for_each_Vertex_around(grid,curVertex2,curParticle,l,m,n){
				secondTerm -= fetch_GradW_withVertex(curParticle,l,m,n)*(fetch_W_withVertex(curParticle,i,j,k)*(curVertex1->A*curVertex2->A));
			}end_for_each_Vertex_around
		}end_for_each_Vertex_around
		curParticle->Momentum += (firstTerm+secondTerm)*deltaT;
	}end_for_each_Particle(curParticle)
}

void CSPIC::updateY(const Range& range){

	static double constAcceletor=deltaT * Square(LIGHT_SPEED)/ (4 * M_PI);

	Vector3D firstTerm,secondTerm;
	Vector3D dFirstTerm,dSecondTerm;
	Vector3D dY,CTCA;
	Vertex *curVertex1;
	Vertex *curVertex2;
	Particle *curParticle;

	int i,j,k;
	int l,m,n;

	for_each_Vertex_within(grid, curVertex1, range){

		firstTerm.restore();
		secondTerm.restore();

		Cache_for_each_Particle_around(grid, curParticle, curVertex1,i,j,k){

			dFirstTerm=curParticle->Momentum;
			dFirstTerm*=fetch_W_withParticle(curParticle,i,j,k);
			firstTerm-=dFirstTerm;

			//firstTerm -= curParticle->Momentum * _W(ParticleRealPosition - curVertex1->_r);

			Cache_for_each_Vertex_around(grid, curVertex2, curParticle,l,m,n){

				dSecondTerm=curVertex2->A;
				dSecondTerm*=fetch_W_withParticle(curParticle,i,j,k)* fetch_W_withVertex(curParticle,l,m,n);
				secondTerm +=  dSecondTerm;

				//cout << "#  " << endl;
				//secondTerm += curVertex2->A * (_W(ParticleRealPosition - curVertex1->_r) * _W(curParticle->Position - Vertex2RealPosition));

			}end_for_each_Vertex_around

		}end_for_each_Particle_around(curParticle)

		CTCA=*curldTCurldA[curVertex1->_x][curVertex1->_y][curVertex1->_z];CTCA*=constAcceletor;

		dY=firstTerm;dY-=secondTerm;
		dY*=deltaT;dY-=CTCA;

		curVertex1->Y+=dY;

		//curVertex1->Y += (firstTerm - secondTerm) * deltaT - *curldTCurldA[curVertex1->_x][curVertex1->_y][curVertex1->_z] * constAcceletor;

	}end_for_each_Vertex_within

}



void CSPIC::updateX(const Range& range){
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

void CSPIC::updateA(const Range& range){
	static double constAcceletor=deltaT * (FOUR_PI);
	Vector3D dA;
	Vertex *curVertex;
	for_each_Vertex_within(grid, curVertex, range){

		dA=curVertex->Y;dA*=constAcceletor;

		curVertex->A += dA;

		//curVertex->A += curVertex->Y * (deltaT * (FOUR_PI));
	}end_for_each_Vertex_within
}

#else

void CSPIC::buidCoefficientTensor(Tensor3D& coefficientTensor,const Range &range, double deltaT, Particle *curParticle){
	coefficientTensor.restore();
	Vertex* curVertex;
	Vector3D colVector;
	Vector3D rowVector;
	Vector3D VertexRealPosition;
	Vector3D acceletor_D;
	for_each_Vertex_around(grid, curVertex, curParticle,VertexRealPosition){

#if OPTIMIZE
		acceletor_D=curParticle->Position;
		acceletor_D-=VertexRealPosition;
		GradW(colVector,acceletor_D);

		rowVector = curVertex->A;
		rowVector *= deltaT;
#else

		GradW(colVector,curParticle->Position - VertexRealPosition);
		rowVector = curVertex->A*deltaT;
#endif


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

void CSPIC::buildRHSVector(Vector3D &RHS,const Range &range, double deltaT, Particle *curParticle){
	RHS.restore();
	Vertex *curVertex1;
	Vertex *curVertex2;
	Vector3D Vertex1RealPosition,Vertex2RealPosition;
	Vector3D acceletor_GradW;
	Vector3D acceletor_D1,acceletor_D2;
	//TODO: can optimize
	for_each_Vertex_around(grid, curVertex1, curParticle,Vertex1RealPosition){
		for_each_Vertex_around(grid, curVertex2, curParticle,Vertex2RealPosition){

#if OPTIMIZE

			acceletor_D1=acceletor_D2=curParticle->Position;
			acceletor_D1-=Vertex1RealPosition;
			acceletor_D2-=Vertex2RealPosition;

			GradW(acceletor_GradW, acceletor_D2);
			acceletor_GradW *= _W(acceletor_D1) * (curVertex1->A * curVertex2->A);
			RHS +=  acceletor_GradW;

#else
			GradW(acceletor_GradW, curParticle->Position-Vertex2RealPosition);
			acceletor_GradW *= _W(curParticle->Position - Vertex1RealPosition) * (curVertex1->A * curVertex2->A);
			RHS +=  acceletor_GradW;
#endif

			//RHS +=GradW(curParticle->Position - Vertex2RealPosition) * (_W(curParticle->Position - Vertex1RealPosition) * (curVertex1->A * curVertex2->A));
		}end_for_each_Vertex_around
	}end_for_each_Vertex_around
	RHS *= -deltaT;
	RHS +=curParticle->Momentum;
}

void CSPIC::updateP(const Range &range){
	Particle *curParticle;

	Tensor3D coefficient;    // The 3*3 coefficient matrix
	Vector3D RHS;                // Right hand side

	for_each_Particle_within(grid, curParticle, range){
		buidCoefficientTensor(coefficient, range, deltaT, curParticle);
		buildRHSVector(RHS, range, deltaT, curParticle);
		rootOfLinearEquationSet(curParticle->Momentum, coefficient, RHS);
	}end_for_each_Particle(curParticle)

}

void CSPIC::updateY(const Range &range){

	Vector3D firstTerm;
	Vector3D secondTerm;
	Vertex *curVertex1;
	Vertex *curVertex2;
	Particle *curParticle;
	Vector3D Vertex2RealPosition,ParticleRealPosition;
	Vector3D acceletor_D1,acceletor_D2,acceletor_D3;
	Vector3D acceletor_FT,acceletor_A;

	static double constAcceletor=deltaT * Square(LIGHT_SPEED)/ (4 * M_PI);

	setCurldA(range);
	setCurldTCurldA(range);

	for_each_Vertex_within(grid, curVertex1, range){
		firstTerm.restore();
		secondTerm.restore();
		for_each_Particle_around(grid, curParticle, curVertex1,ParticleRealPosition){

#if OPTIMIZE
			acceletor_D1=ParticleRealPosition;
			acceletor_D1-=curVertex1->_r;

			acceletor_FT=curParticle->Momentum;
			acceletor_FT*=_W(acceletor_D1);

			firstTerm-=acceletor_FT;

#else
			firstTerm -= curParticle->Momentum * _W(ParticleRealPosition - curVertex1->_r);
#endif

			/////Vertex2RealPosition is not real here.
			for_each_Vertex_around(grid, curVertex2, curParticle,Vertex2RealPosition){

#if OPTIMIZE
				acceletor_D2=ParticleRealPosition;
				acceletor_D2-=curVertex1->_r;

				//as we only use relative Position, here we use original Particle position to recorrect the 'real' vertex position.
				acceletor_D3=curParticle->Position;
				acceletor_D3-=Vertex2RealPosition;

				acceletor_A=curVertex2->A;
				acceletor_A*=_W(acceletor_D2) * _W(acceletor_D3);

				secondTerm +=  acceletor_A;

#else
				secondTerm += curVertex2->A * (_W(ParticleRealPosition - curVertex1->_r) * _W(curParticle->Position - Vertex2RealPosition));
#endif

			}end_for_each_Vertex_around
		}end_for_each_Particle_around(curParticle)
		curVertex1->Y += (firstTerm - secondTerm) * deltaT - *curldTCurldA[curVertex1->_x][curVertex1->_y][curVertex1->_z] * constAcceletor;

	}end_for_each_Vertex_within

}

void CSPIC::updateX(const Range &range){
	Particle *curParticle;

	Vertex *curVertex;
	Vector3D secondTerm;
	Vector3D accelerator_A;
	Vector3D accelerator_P;
	Vector3D VertexRealPosition;
	Vector3D accelerator_D;

	for_each_Particle_within(grid, curParticle, range){
		secondTerm.restore();

		for_each_Vertex_around(grid, curVertex,curParticle, VertexRealPosition){
#if OPTIMIZE

			accelerator_D=curParticle->Position;
			accelerator_D-=VertexRealPosition;

			accelerator_A=curVertex->A;
			accelerator_A*=_W(accelerator_D);

			secondTerm += accelerator_A;

#else
			secondTerm += curVertex->A * _W(curParticle->Position - VertexRealPosition);
#endif
		}end_for_each_Vertex_around

#if OPTIMIZE
		accelerator_P=curParticle->Momentum;
		accelerator_P+=secondTerm;
		accelerator_P*=deltaT;

		curParticle->Position += accelerator_P;
#else
		curParticle->Position += (curParticle->Momentum + secondTerm) * deltaT;
#endif

	}end_for_each_Particle(curParticle)

}

void CSPIC::updateA(const Range &range){
	Vertex *curVertex;
	for_each_Vertex_within(grid, curVertex, range){
		curVertex->A += curVertex->Y * (deltaT * (FOUR_PI));

	}end_for_each_Vertex_within
}
#endif
