/*
 * Distributer.cpp
 *
 *  Created on: Apr 6, 2016
 *      Author: zlstudio
 */

#include <Buffer.h>
#include "Distributer.h"

Distributer::Distributer(int _routes,int _thisroute,Grid* mananged_grid){
	routes=_routes;
	thisRoute=_thisroute;

	grid=mananged_grid;

	sendBuffer=new double*[_routes];

	sendLength=new int[_routes];

	particleBuffer=new Cell[_routes]();

}

Distributer::~Distributer() {
}

void Distributer::pack(Particle *p,int dest){

	particleBuffer[dest].pushParticle(p);

}

void Distributer::assembleAndDispatch(){
	Particle* p;
	for (int r = 0; r < routes; ++r) {

		sendLength[r]=particleBuffer[r].length();

		if(r==thisRoute||particleBuffer[r].length()==0)continue;

		sendBuffer[r]=new double[sendLength[r]*6];

		buffer_init();

		for_each_Particle_in_Cell((particleBuffer+r),p){
			buffer_Write_Vector3D(sendBuffer[r],p->Position);
			buffer_Write_Vector3D(sendBuffer[r],p->Momentum);
		}end_for_each_Particle_in_Cell(p)

//cout<<"send:"<<sendBuffer[r][0]<<","<<sendBuffer[r][1]<<","<<sendBuffer[r][2]<<endl;
		particleBuffer[r].destroyAllParticlesInside();

	}

	//Send & Recv in stair-like process
	/**S=Send; R=Recv
	 * 1:          \ S-2,3,4 \ R-2,3,4
	 * 2: R-1     \ S-1,3,4 \ R-3,4
	 * 3: R-2,3   \ S-1,2,4 \ R-4
	 * 4: R-2,3,4 \ S-1,2,3 \
	 */

	//	trace(thisRoute<<" waiting 0~"<<thisRoute)

	for (int i = 0; i < thisRoute; ++i) {
		collect(i);
	}

	for (int i = 0; i < routes; ++i) {
		if (i==thisRoute)continue;
		dispatch(i);
	}

	for (int i = thisRoute+1; i < routes; ++i) {
		collect(i);
	}

}


void Distributer::dispatch(int dest){

	send_signal_to(dest,sendLength[dest]);

	if(sendLength[dest]==0)return;

	//trace(thisRoute<<" send "<<sendLength[dest]<<" to " <<dest);

	buffer_send_to(sendBuffer[dest],sendLength[dest]*6,dest,SYNC_PARTICLE);

	delete sendBuffer[dest];

}

void Distributer::collect(int dest){

	if (dest==thisRoute)return;

	 Electron* e_in;

	//get the length of comming particles
	int incoming=wait_signal_from(dest);

	if(incoming==0)return;

	//receive data
	buffer_fast_init(recvbuffer,incoming*6);
	buffer_fast_recv_from(recvbuffer,dest,SYNC_PARTICLE);

	for (int j = 0; j < incoming; ++j) {
		e_in=new Electron();
		buffer_fast_Read_Vector3D(recvbuffer,e_in->Position);
		buffer_fast_Read_Vector3D(recvbuffer,e_in->Momentum);

		grid->addParticle(e_in);
	}

}


