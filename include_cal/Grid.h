/*
 * Grid.h
 *
 *  Created on: Mar 3, 2016
 *      Author: zlstudio
 */


#ifndef GRID_H_
#define GRID_H_

#include <Range.h>
#include <Vertex.h>
#include <Cell.h>
#include <PhysicalParameters.h>
#include "Particle.h"
#include "Range.h"


#define trace(inf) cout<<inf<<endl;

#define for_each_Cell(grid1,cell_iterator) \
		for (int _i = 0; _i < grid1->gridX(); ++_i) {\
			for (int _j = 0; _j < grid1->gridY(); ++_j) {\
				for (int _k = 0; _k < grid1->gridZ(); ++_k) {\
					cell_iterator=grid1->cell(_i,_j,_k);

#define end_for_each_Cell }}}

#define for_each_Cell_within(grid1,iterator,range) \
		for (int _i = range.sx; _i < range.ex; ++_i) {\
			for (int _j = range.sy; _j < range.ey; ++_j) {\
				for (int _k = range.sz; _k < range.ez; ++_k) {\
					iterator=grid1->cell(_i,_j,_k);

#define for_each_Vertex(grid1,iterator) \
		for (int _i = 0; _i < grid1->gridX(); ++_i) {\
			for (int _j = 0; _j < grid1->gridY(); ++_j) {\
				for (int _k = 0; _k < grid1->gridZ(); ++_k) {\
					iterator=grid1->vertex(_i,_j,_k);

#define end_for_each_Vertex }}}

#define for_each_Vertex_within(grid1,iterator,range) \
		for (int _i = range.sx; _i < range.ex; ++_i) {\
			for (int _j = range.sy; _j < range.ey; ++_j) {\
				for (int _k = range.sz; _k < range.ez; ++_k) {\
					iterator=grid1->vertex(_i,_j,_k);

#define end_for_each_Vertex_within }}}

#if OPTIMIZE_1D_MPI
#define for_each_Vertex_around(grid1,iterator,particle,realposition) \
		int sx=particle->cell->_x-1,sy=particle->cell->_y-1,sz=particle->cell->_z-1;\
		int ex=particle->cell->_x+2,ey=particle->cell->_y+2,ez=particle->cell->_z+2;\
		for (int _l = sx; _l <= ex; _l++){\
			for (int _m = sy; _m <=ey ;_m++){\
				for (int _n = sz; _n <= ez; _n++){\
					iterator=grid->vertex(_l, _m, _n);\
					realposition.x=iterator->_x;\
					realposition.y=iterator->_y+_m;\
					realposition.z=iterator->_z+_n;

#elif OPTIMIZE_1D
//As getPeriod(i,1)==>i, gridY()=gridZ()==>1. gridX()>=2, then getPeriod will be confine to -1,0,1
#define for_each_Vertex_around(grid1,iterator,particle,realposition) \
		int sx=particle->cell->_x-1,sy=particle->cell->_y-1,sz=particle->cell->_z-1;\
		int ex=particle->cell->_x+2,ey=particle->cell->_y+2,ez=particle->cell->_z+2;\
		for (int _l = sx; _l <= ex; _l++){\
			for (int _m = sy; _m <=ey ;_m++){\
				for (int _n = sz; _n<= ez; _n++){\
					iterator=grid->vertex(_l, _m, _n);\
					realposition.x=iterator->_x+Fast_adjustLoop(_l,grid->_width);\
					realposition.y=iterator->_y+_m;\
					realposition.z=iterator->_z+_n;

#elif OPTIMIZE_3D
//when x&y&z are all greater than 1, then getPeriod will be confine to -1,0,1
#define for_each_Vertex_around(grid1,iterator,particle,realposition) \
		int sx=particle->cell->_x-1,sy=particle->cell->_y-1,sz=particle->cell->_z-1;\
		int ex=particle->cell->_x+2,ey=particle->cell->_y+2,ez=particle->cell->_z+2;\
		for (int _l = sx; _l <= ex; _l++){\
			for (int _m = sy; _m <=ey ;_m++){\
				for (int _n = sz; _n<= ez; _n++){\
					iterator=grid->vertex(_l, _m, _n);\
					realposition.x=iterator->_x+Fast_adjustLoop(_l,grid->_width);\
					realposition.y=iterator->_y+Fast_adjustLoop(_m,grid->_height);\
					realposition.z=iterator->_z+Fast_adjustLoop(_n,grid->_length);

#else
#define for_each_Vertex_around(grid1,iterator,particle,realposition) \
		int sx=particle->cell->_x-1,sy=particle->cell->_y-1,sz=particle->cell->_z-1;\
		int ex=particle->cell->_x+2,ey=particle->cell->_y+2,ez=particle->cell->_z+2;\
		for (int _l = sx; _l <= ex; _l++){\
			for (int _m = sy; _m <=ey ;_m++){\
				for (int _n = sz; _n<= ez; _n++){\
					iterator=grid->vertex(_l, _m, _n);\
					realposition.x=iterator->_x+getPeriod(_l,grid->gridX())*grid->gridX();\
					realposition.y=	iterator->_y+getPeriod(_m,grid->gridY())*grid->gridY();\
					realposition.z=	iterator->_z+getPeriod(_n,grid->gridZ())*grid->gridZ();
#endif

#define end_for_each_Vertex_around }}}


#if OPTIMIZE_1D_MPI
#define for_each_Particle_around(grid1,iterator,vertex,realposition) \
		int sx=vertex->_x-2,sy=vertex->_y-2,sz=vertex->_z-2;\
		int ex=vertex->_x+1,ey=vertex->_y+1,ez=vertex->_z+1;\
		for (int _l = sx; _l <= ex; _l++){\
			for (int _m = sy; _m <=ey;_m++){\
				for (int _n = sz; _n<= ez; _n++){\
					iterator=grid->cell(_l, _m, _n)->firstParticle();\
					while(iterator!=NULL){\
						realposition.x=iterator->Position.x;\
						realposition.y=iterator->Position.y+_m;\
						realposition.z=iterator->Position.z+_n;

#elif OPTIMIZE_1D
//As getPeriod(i,1)==>i, here gridY()=gridZ()==>1, thus gerPeriod(i,gridY())*gridY() = i. gridX()>=2, then getPeriod will be confine to -1,0,1
#define for_each_Particle_around(grid1,iterator,vertex,realposition) \
		int sx=vertex->_x-2,sy=vertex->_y-2,sz=vertex->_z-2;\
		int ex=vertex->_x+1,ey=vertex->_y+1,ez=vertex->_z+1;\
		for (int _l = sx; _l <= ex; _l++){\
			for (int _m = sy; _m <=ey;_m++){\
				for (int _n = sz; _n<= ez; _n++){\
					iterator=grid->cell(_l, _m, _n)->firstParticle();\
					while(iterator!=NULL){\
						realposition.x=iterator->Position.x+Fast_adjustLoop(_l,grid->_width);\
						realposition.y=iterator->Position.y+_m;\
						realposition.z=iterator->Position.z+_n;

#elif OPTIMIZE_3D
//when x&y&z are all greater than 1, then getPeriod will be confine to -1,0,1
#define for_each_Particle_around(grid1,iterator,vertex,realposition) \
		int sx=vertex->_x-2,sy=vertex->_y-2,sz=vertex->_z-2;\
		int ex=vertex->_x+1,ey=vertex->_y+1,ez=vertex->_z+1;\
		for (int _l = sx; _l <= ex; _l++){\
			for (int _m = sy; _m <=ey;_m++){\
				for (int _n = sz; _n<= ez; _n++){\
					iterator=grid->cell(_l, _m, _n)->firstParticle();\
					while(iterator!=NULL){\
						realposition.x=iterator->Position.x+Fast_adjustLoop(_l,grid->_width);\
						realposition.y=iterator->Position.y+Fast_adjustLoop(_m,grid->_height);\
						realposition.z=iterator->Position.z+Fast_adjustLoop(_n,grid->_length);

#else
#define for_each_Particle_around(grid1,iterator,vertex,realposition) \
		int sx=vertex->_x-2,sy=vertex->_y-2,sz=vertex->_z-2;\
		int ex=vertex->_x+1,ey=vertex->_y+1,ez=vertex->_z+1;\
		for (int _l = sx; _l <= ex; _l++){\
			for (int _m = sy; _m <=ey;_m++){\
				for (int _n = sz; _n<= ez; _n++){\
					iterator=grid->cell(_l, _m, _n)->firstParticle();\
					while(iterator!=NULL){\
						realposition.x=iterator->Position.x+getPeriod(_l,grid->gridX())*grid->gridX();\
						realposition.y=	iterator->Position.y+getPeriod(_m,grid->gridY())*grid->gridY();\
						realposition.z=	iterator->Position.z+getPeriod(_n,grid->gridZ())*grid->gridZ();

#endif

#define end_for_each_Particle_around(iterator) iterator=iterator->nextParticle;}}}}

#define for_each_Particle(grid1,iterator) \
		for (int _i = 0; _i < grid1->gridX(); ++_i) {\
			for (int _j = 0; _j < grid1->gridY(); ++_j) {\
				for (int _k = 0; _k < grid1->gridZ(); ++_k) {\
					iterator=grid1->cell(_i,_j,_k)->firstParticle();\
					while(iterator!=NULL){

#define end_for_each_Particle(iterator) iterator=iterator->nextParticle;}}}}

#define end_for_each_Particle_raw }}}}


///////////
#define for_each_Particle_within(grid1,iterator,range) \
		for (int _i = range.sx; _i < range.ex; ++_i) {\
			for (int _j = range.sy; _j < range.ey; ++_j) {\
				for (int _k = range.sz; _k < range.ez; ++_k) {\
					iterator=grid1->cell(_i,_j,_k)->firstParticle();\
					while(iterator!=NULL){

///////////
#define in_Range(i,j,k,range) \
		for (i = range.sx; i < range.ex; ++i) {\
			for (j = range.sy; j < range.ey; ++j) {\
				for (k = range.sz; k < range.ez; ++k) {\

#define end_in_Range }}}


////////Cal Period
#define getPeriod(i,loop) ((i>=0)?(i)/loop:(i-loop+1)/loop)

//when Period multiplies loop
#define adjustLoop(i,loop) ((i>=0)?(i)/loop*loop:(i-loop+1)/loop*loop)

//for circumstance confined to period -1,0,1 and then multiplies loop. when loop==1, this should be replaced by i.
#define Fast_adjustLoop(i,loop) ((i<0)?(-loop):((i>=loop)?loop:0))

//(i+loop)%loop==> i<0 i+=loop; i>=loop i-=loop
//to get position in the period one locates. when period is confined to -1,0,1
#define Fast_adjustPosition(i,loop) ((i<0)?(i+loop):((i>=loop)?i-loop:i))


/*
 * Grid=Vertex+Cell
 */

class Vertex;
class Cell;
class Range;

class Grid {
	public:
		Grid(int x, int y,int z, bool period);
		virtual ~Grid();

		Cell* cell(int a,int b ,int c);
		Vertex* vertex(int a, int b, int c);


		////KEY FUNC

		bool addParticle(Particle* p);

		virtual bool directAddParticle(Particle* p);

		bool moveParticleTo(Particle* p, const Vector3D &location);

		virtual void refreshParticleLocation();

		virtual void syncSharedParticles();

		virtual void reportToHost();

		////

		int gridX();
		int gridY();
		int gridZ();

		int particles();

		void showGridMap();

		Vector3D**** A_indexer;
		Vector3D**** Y_indexer;


		///obsolete
		int lengthX();
		int lengthY();
		int lengthZ();

		double scale();

		//sacrifies

		Range World;
		Range workSpace;
		int _width,_height,_length;
		int _workLength;

		int workLength();

	protected:
		void initGrid3D(int x,int y,int z);

		Vertex**** VertexContainer;
		Cell**** CellContainer;

		Cell* particleSwap;

		double _scale;

		int _particle_amount;

		bool _period;
};

#endif /* GRID_H_ */
