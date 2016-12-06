/*
 * indexCache.h
 *
 *  Created on: 2016-5-10
 *      Author: zlstudio
 */

#ifndef INCLUDE_INDEXCACHE_H_
#define INCLUDE_INDEXCACHE_H_


#if USE_CACHE


///////Vertex Around
//Particles need vertexes around its cell:
//-1,0,p,1,2
//to use index, we have to +1.

#include <Particle.h>
#include <Vertex.h>
#include <Grid.h>
#include <Vector3D.h>
#include <CSPIC.h>



///use cache
///we need the vertex so we have to keep it, and we now need index 0~3 instead of realPosition
#define Cache_for_each_Vertex_around(grid1,iterator,particle,ll,mm,nn) \
		int sx=particle->cell->x()-1,sy=particle->cell->y()-1,sz=particle->cell->z()-1;\
		for (ll = 0; ll <= 3; ll++){\
			for (mm = 0; mm <=3 ;mm++){\
				for (nn = 0; nn <= 3; nn++){\
					iterator=grid->vertex(ll+sx, mm+sy, nn+sz);\

#define fetch_W_withVertex(p,ll,mm,nn) p->W_cache[ll][mm][nn]
#define fetch_GradW_withVertex(p,ll,mm,nn) p->GW_cache[ll][mm][nn]





///////////Particle Around



//#define Cache_for_each_Particle_around(grid1,iterator,vertex,ll,mm,nn)\
//		int sx=vertex->_x-2,sy=vertex->_y-2,sz=vertex->_z-2;\
//		for (ll = 0; ll <= 3; ll++){\
//			for (mm = 0; mm <=3;mm++){\
//				for (nn = 0; nn <= 3; nn++){\
//					iterator=grid->cell(ll+sx, mm+sy, nn+sz)->firstParticle();\
//					while(iterator!=NULL){

#define Cache_for_each_Particle_around(grid1,iterator,vertex,ll,mm,nn)\
		int sx=vertex->_x-2,sy=vertex->_y-2,sz=vertex->_z-2;\
		for (ll = 0; ll <= 3; ll++){\
			for (mm = 2; mm <=2;mm++){\
				for (nn = 2; nn <= 2; nn++){\
					iterator=grid->cell(ll+sx, mm+sy, nn+sz)->firstParticle();\
					while(iterator!=NULL){

#define fetch_W_withParticle(p,ll,mm,nn) p->W_cache[3-ll][3-mm][3-nn]
#define fetch_GradW_withParticle(p,ll,mm,nn) p->GW_cache[3-ll][3-mm][3-nn]




#endif


#endif /* INCLUDE_INDEXCACHE_H_ */
