/*
 * Eular.cpp
 *
 *  Created on: Mar 6, 2016
 *      Author: zlstudio
 */

#include <Eular.h>

Grid* Eular::grid;

inline void delayMS(int length){
	struct timespec sleeper, dummy ;
	sleeper.tv_sec  = (time_t)(length / 1000) ;
	sleeper.tv_nsec = (long)(length % 1000) * 1000000 ;
	nanosleep (&sleeper, &dummy) ;
}

void Eular::nextStep(const Range &range,double dt){
	Particle*  particle_iterator;
	for_each_Particle_within(grid,particle_iterator,range){
		particle_iterator->Position+=particle_iterator->Momentum*dt;
	}end_for_each_Particle(particle_iterator);

	//grid->showGridMap();
}

void Eular::nextStep2(const Range &range,double dt){
	Particle *p1,*p2;

	for_each_Particle(grid,p1){

		Vector3D speed(0,0,0);
		double d;

		for_each_Particle(grid,p2){

			if(p2!=p1){
				d=p2->Position>>p1->Position;
				speed+= (p1->Position-p2->Position)/(d*d*d/p1->charge/p2->charge);
			}

		}end_for_each_Particle(p2);

		speed*=dt/p1->mass;

		p1->Momentum+=speed;

	}end_for_each_Particle(p1);


	for_each_Particle(grid,p1){

		grid->moveParticleTo(p1,p1->Position+p1->Momentum*dt);

	}end_for_each_Particle(p1);

}

void Eular::nextStep3(const Range &range,double dt){
	updateY(range,dt);
	updateP(range,dt);
#pragma omp barrier
	updateX(range,dt);
}

void Eular::updateP(const Range &range,double dt){
	Particle  *p,*center;

	for_each_Particle_within(grid,center,range){
		Vector3D speed(0,0,0);
		double d;

		int range=2;

		int sx=_i-range;
		int sy=_j-range;
		int sz=_k-range;

		int ex=_i+range+1;
		int ey=_j+range+1;
		int ez=_k+range+1;

		Range r(sx,ex,sy,ey,sz,ez);

		for_each_Particle_within(grid,p,r){

			if(p!=center){

				Vector3D realPosition;

				realPosition.x=p->Position.x+getPeriod(_i,grid->gridX())*grid->gridX();
				realPosition.y=p->Position.y+getPeriod(_j,grid->gridY())*grid->gridY();
				realPosition.z=p->Position.z+getPeriod(_k,grid->gridZ())*grid->gridZ();

				d=realPosition>>center->Position;

				speed+= (center->Position-realPosition)/(d*d*d);
			}

		}end_for_each_Particle(p)

		speed*=dt;

		center->Momentum+=speed;

	}end_for_each_Particle(center)
}

void Eular::updateX(const Range &range,double dt){
	Particle  *center;

	for_each_Particle_within(grid,center,range){
		center->Position+=center->Momentum*dt;
	}end_for_each_Particle(center)

}

void Eular::updateY(const Range &range,double dt){
	Particle  *p;
	Vertex *v;
	double d;

	Vector3D realPosition;

	for_each_Vertex_within(grid,v,range){
		Vector3D fieldE(0,0,0);

		for_each_Particle_around(grid,p,v,realPosition){

			d=realPosition>>v->r();

			fieldE+= (p->Position-realPosition)/(d*d*d);

		}end_for_each_Particle_around(p)

		v->Y=fieldE;

	}end_for_each_Vertex

}


void Eular::nextStep4(const Range &range,double dt){
	Particle*  particle_iterator;

	for_each_Particle_within(grid,particle_iterator,range){

		particle_iterator->Momentum+=grid->vertex(_i,_j,_k)->A*dt;

	}end_for_each_Particle(particle_iterator);


}


