/*
 * GridManager.cpp
 *
 *  Created on: Apr 6, 2016
 *      Author: zlstudio
 */

#include <MPIGrid.h>
#include <omp.h>
#include <RunManager.h>
#include "Buffer.h"
#include "Distributer.h"

////////MPI_GRID architecture
////////prevNode <-- thisNode --> nextNode
////////Architecture for each Node:

////////|	   	 	 	 			storeLength 	 			 		 	|///////worklength+6
////////|      SWAP	|				workSpace				|     SWAP	|///////workSpace&SWAP
////////|       VadjL	|        VshL	|				|       VshR	|       VadjR	|///////Vertex shared(sh)/adjacent(adj) space x3
////////|	|   PadjL	|   PshL 	|						|   PshR	|   PadjR	|	|///////Particle shared(sh)/adjacent(adj) space x2
////////| SWAP_LEN	|<==SWAP_LEN	storeLength-SWAP_LEN==>	| SWAP_LEN	|

MPIGrid::MPIGrid(int w,int h,int l):
Grid(w/RunManager::Nodes+SWAP_LEN*2,h,l,true),
nodes(RunManager::Nodes),thisNode(RunManager::MPI_ID),
_worldwidth(w),_worldheight(h),_worldlength(l),
_storeLength(_worldwidth/nodes+SWAP_LEN*2),

VxAdjacentL(0,SWAP_LEN,0,_worldheight,0,_worldlength),
VxSharedL(SWAP_LEN,V_SHARE_END,0,_worldheight,0,_worldlength),
VxSharedR(_storeLength-V_SHARE_END,_storeLength-SWAP_LEN,0,_worldheight,0,_worldlength),
VxAdjacentR(_storeLength-SWAP_LEN,_storeLength,0,_worldheight,0,_worldlength),

PtAdjacentL(P_ADJ_BEGIN_L,SWAP_LEN,0,_worldheight,0,_worldlength),
PtSharedL(SWAP_LEN,P_SHARE_END_L,0,_worldheight,0,_worldlength),
PtSharedR(_storeLength-P_SHARE_END_R,_storeLength-SWAP_LEN,0,_worldheight,0,_worldlength),
PtAdjacentR(_storeLength-SWAP_LEN,_storeLength-P_ADJ_BEGIN_R,0,_worldheight,0,_worldlength)
{

	World=Range(0,w,0,h,0,l);
	workSpace=Range(SWAP_LEN,_storeLength-SWAP_LEN,0,_worldheight,0,_worldlength);
	//Init MPI cal grid
	if(_worldwidth%nodes!=0)cerr<<"Warning: Grid size does not match Network architecture!"<<endl;

	//Prev & Next & This
	prevNode=(thisNode==0?nodes-1:thisNode-1);
	nextNode=(thisNode==nodes-1?0:thisNode+1);

	distributer=new Distributer(RunManager::Nodes,RunManager::MPI_ID,this);

	_workLength=_worldwidth/nodes;
#if REPORT
	report_socket=connectSocket(HOST_ADDRESS,TRANSMITTER_IN_PORT);
	send(report_socket,&RunManager::MPI_ID,sizeof(RunManager::MPI_ID),0);
#endif
}

MPIGrid::~MPIGrid() {
	// TODO Auto-generated destructor stub
}


bool MPIGrid::directAddParticle(Particle* p){

	p->Position.x+=SWAP_LEN;

	if(CellContainer[3][0][0]->pushParticle(p)){
		_particle_amount++;
		return true;
	}

	cerr<<"Warning from MPIGrid : Adding Particle Failed!"<<endl;

	return false;
}

void MPIGrid::refreshParticleLocation(){
	//	static int n=0;
	//static double dt1=0,dt2=0,dt3=0;

	//	double t=omp_get_wtime();
	syncVertex();
	//	dt2+=omp_get_wtime()-t;

	//	t=omp_get_wtime();
	syncParticles();
	//dt1+=omp_get_wtime()-t;

	//	t=omp_get_wtime();
	shareParticles();
	//dt3+=omp_get_wtime()-t;

	//n++;

	//	if(thisNode==0&&n>=10){
	//		cout<<"SYNC Time:\tVertex\t"<<dt2*100<<" ms\tParticles:\t"<<dt1*100<<" ms\t Share:\t"<<dt3*100<<endl;
	//		dt1=dt2=dt3=n=0;
	//	}
}

void MPIGrid::report(){
	//TODO: report to node 0
}

void MPIGrid::reportToHost(){


	reportParticles();
	reportVertexes();

}
//void MPIGrid::addParticle(Particle* p){
//	//Trans-coordinate to local
//	p->Position.x+=2;
//	p->Position.x-=thisNode*workLength;
//
//	Grid::addParticle(p);
//}

void MPIGrid::reportParticles(){
#if REPORT

	int n=0;

	Particle* p;

	Cell* c;

	for_each_Cell_within(this,c,workSpace){

		n+=c->length();

	}end_for_each_Cell;

	buffer_fast_init(buffer,n*6);

	for_each_Particle_within(this,p,workSpace){

		buffer_fast_Write_Vector3D(buffer,p->Position);
		buffer_fast_Write_Vector3D(buffer,p->Momentum);

	}end_for_each_Particle(p);

	send(report_socket,&n,sizeof(n),0);

	send(report_socket,buffer,sizeof(buffer),0);
#endif
}

void MPIGrid::reportVertexes(){
#if REPORT

	Vertex* v;

	int V=workSpace.Volumn();

	buffer_fast_init(buffer,V*6);

	for_each_Vertex_within(this,v,workSpace){

		buffer_fast_Write_Vector3D(buffer,v->A);
		buffer_fast_Write_Vector3D(buffer,v->Y);

	}end_for_each_Vertex

	send(report_socket,&V,sizeof(V),0);

	send(report_socket,buffer,sizeof(buffer),0);
#endif
}


void MPIGrid::syncVertex(){

	/////shared Right====>external Left
	if(!(thisNode%2))sendVertexes(nextNode,VxSharedR);

	recvVertexes(prevNode,VxAdjacentL);

	if(thisNode%2)sendVertexes(nextNode,VxSharedR);

	/////shared Left====>external Right
	if(!(thisNode%2))sendVertexes(prevNode,VxSharedL);

	recvVertexes(nextNode,VxAdjacentR);

	if(thisNode%2)sendVertexes(prevNode,VxSharedL);

}

void MPIGrid::syncParticles(){

	Particle* p,*p_swap;

	//at this point of time, all particles previously locate in work space still exist in the list of original cell.
	//So we need only iterate within workSpace

	for_each_Particle_within(this,p,workSpace){

		if(p->cell!=NULL){

			//get the real local X
			double real_x=p->Position.x-SWAP_LEN;

			//we have to use floor as double -> int make negative number bigger.
			int real_cell_x=floor(real_x);

			//We can do int=double because after switch period these will certainly greater than 0;
			real_cell_x=real_x-=getPeriod(real_cell_x,_worldwidth)*_worldwidth;

			p->Position.x=real_x+SWAP_LEN;

			///local cell
			int cell_x=real_cell_x+SWAP_LEN;
			int cell_y=floor(p->Position.y);
			int cell_z=floor(p->Position.z);

			//We can do int=double because after switch period these will certainly greater than 0;
			cell_y=p->Position.y-=getPeriod(cell_y,_worldheight)*_worldheight;
			cell_z=p->Position.z-=getPeriod(cell_z,_worldlength)*_worldlength;


			if(real_cell_x>=_workLength){ //here real_cell_x would not be negative. 0 is the begin of the local grid.
				///Out of boundary --> to be sent

				p_swap=p;

				p=p->nextParticle;

				//this Node is alway 0 to the local while real_cell_x is always positive, so we can get relate Node as below:
				int relateNodePath=real_cell_x/_workLength;

				//Trans to destination coordinate
				p_swap->Position.x-=relateNodePath*_workLength;

				//destination:
				int dest=thisNode+relateNodePath;

				//if it is 3 -> 4 and nodes=4, we will have to cast 3-nodes. dest is surely already positive.
				//of course, here dest has already be restrict in +- 1 period.
				dest-=(dest>=nodes)*nodes;

				p_swap->cell->withdrawParticle(p_swap);
				--_particle_amount;

				//pack up
				distributer->pack(p_swap,dest);

			}else if(cell(cell_x,cell_y,cell_z)!=p->cell){

				///Inside this Grid, but in a different cell

				p_swap=p;

				p=p->nextParticle;

				p_swap->cell->withdrawParticle(p_swap);

				//we have to push this particle into Swap first in order to complete the iteration, or some of the particles will be ignored.

				particleSwap->pushParticle(p_swap);

			}else{
				///then the particle need not changing cell
				p=p->nextParticle;
			}
		}else{
			cerr<<"Warning from Grid : Particle with no cell pointer exists!"<<endl;
		}

	}end_for_each_Particle_raw;

	//put the swap particles to where they should go.

	for_each_Particle_in_Cell(particleSwap,p){

		p_swap=p;

		p=p->nextParticle;

		particleSwap->withdrawParticle(p_swap);
		cell(p_swap->Position.x,p_swap->Position.y,p_swap->Position.z)->pushParticle(p_swap);

	}end_for_each_Particle_in_Cell_raw;

	////Particals will be packed for each Node.
	////X & P are going to be packed, 6 doubles for each particle
	distributer->assembleAndDispatch();

}

void MPIGrid::shareParticles(){
	/////shared Right====>external Left

	if(!(thisNode%2)){sendParticles(nextNode,PtSharedR);}

	recvParticles(prevNode,PtAdjacentL);

	if(thisNode%2){sendParticles(nextNode,PtSharedR);}

	/////shared Left====>external Right
	if(!(thisNode%2)){sendParticles(prevNode,PtSharedL);}

	recvParticles(nextNode,PtAdjacentR);

	if(thisNode%2){sendParticles(prevNode,PtSharedL);}

	//trace("Node "<<thisNode<<" shared");
}

void MPIGrid::syncSharedParticles(){

	///for only sync P in adj area, after updateP
	/////shared Right====>external Left

	if(!(thisNode%2)){sendParticles(nextNode,PtSharedR);}

	recvMomentum(prevNode,PtAdjacentL);

	if(thisNode%2){sendParticles(nextNode,PtSharedR);}

	/////shared Left====>external Right
	if(!(thisNode%2)){sendParticles(prevNode,PtSharedL);}

	recvMomentum(nextNode,PtAdjacentR);

	if(thisNode%2){sendParticles(prevNode,PtSharedL);}

	//trace("Node "<<thisNode<<" shared");
}

void MPIGrid::sendParticles(int dest,const Range& r){

	int n=0;

	Particle* p;

	Cell* c;

	for_each_Cell_within(this,c,r){

		n+=c->length();

	}end_for_each_Cell;

	buffer_fast_init(buffer,n*6);

	for_each_Particle_within(this,p,r){

		buffer_fast_Write_Vector3D(buffer,p->Position);
		buffer_fast_Write_Vector3D(buffer,p->Momentum);

	}end_for_each_Particle(p);

	send_signal_to(dest,n);

	buffer_fast_send_to(buffer,dest,SHARE_PARTICLE);
	//MPItrace("send "<<n<<" to "<<dest);

}

void MPIGrid::recvParticles(int source,const Range& r){
	//about 80~100ms ,Node0 20ms
	Cell* c;

	for_each_Cell_within(this,c,r){

		_particle_amount-=c->length();

		c->destroyAllParticlesInside();

	}end_for_each_Cell;

	//then add;

	Electron* e_in;

	int n=wait_signal_from(source);

	//MPItrace("get "<<n);

	buffer_fast_init(recvbuffer,n*6);

	buffer_fast_recv_from(recvbuffer,source,SHARE_PARTICLE);

	for (int i = 0; i < n; ++i) {

		e_in=new Electron();
		buffer_fast_Read_Vector3D(recvbuffer,e_in->Position);
		buffer_fast_Read_Vector3D(recvbuffer,e_in->Momentum);

		//modify the position to local coordinate here.
		source==prevNode?e_in->Position.x-=_workLength:e_in->Position.x+=_workLength;

		addParticle(e_in);
	}
}


void MPIGrid::recvMomentum(int source,const Range& r){

	Particle* p;
	Vector3D tmp;

	int n=wait_signal_from(source);

	//MPItrace("get "<<n);

	buffer_fast_init(recvbuffer,n*6);

	buffer_fast_recv_from(recvbuffer,source,SHARE_PARTICLE);

	for_each_Particle_within(this,p,r){

		buffer_fast_Read_Vector3D(recvbuffer,tmp);
		if(tmp.y!=p->Position.y)cerr<<"Particle does not match!"<<endl;
		buffer_fast_Read_Vector3D(recvbuffer,p->Momentum);


	}end_for_each_Particle(p)

}
//void MPIGrid::recvParticles(int source,const Range& r){
////about the same
//	//We must make use of current existed particles.
//	Cell* c;
//
//	for_each_Cell_within(this,c,r){
//
//		c->restore();
//
//	}end_for_each_Cell;
//
//	//then add;
//
//	Electron p_receive;
//
//	int n=wait_signal_from(source);
//
//	buffer_fast_init(recvbuffer,n*6);
//
//	buffer_fast_recv_from(recvbuffer,source,SHARE_PARTICLE);
//
//	int cell_x,cell_y,cell_z;
//
//	for (int i = 0; i < n; ++i) {
//
//		buffer_fast_Read_Vector3D(recvbuffer,p_receive.Position);
//		buffer_fast_Read_Vector3D(recvbuffer,p_receive.Momentum);
//
//		//modify the position to local coordinate here.
//		source==prevNode?p_receive.Position.x-=workLength:p_receive.Position.x+=workLength;
//
//		cell_x=p_receive.Position.x;
//		cell_y=p_receive.Position.y;
//		cell_z=p_receive.Position.z;
//
//		_particle_amount+=CellContainer[cell_x][cell_y][cell_z]->overwriteParticle(&p_receive);
//	}
//
//	for_each_Cell_within(this,c,r){
//
//		_particle_amount-=c->cut();
//
//	}end_for_each_Cell;
//}

void MPIGrid::sendVertexes(int dest,const Range& r){

	Vertex* vertex;

	int V=r.Volumn()*6;

	buffer_fast_init(sendBuffer,V);

	for_each_Vertex_within(this,vertex,r){

		buffer_fast_Write_Vector3D(sendBuffer,vertex->A);
		buffer_fast_Write_Vector3D(sendBuffer,vertex->Y);

	}end_for_each_Vertex_within

	buffer_fast_send_to(sendBuffer,dest,SYNC_VERTEX);

}

void MPIGrid::recvVertexes(int source,const Range& r){
	Vertex* vertex;

	int V=r.Volumn()*6;

	buffer_fast_init(recvBuffer,V);

	buffer_fast_recv_from(recvBuffer,source,SYNC_VERTEX);

	for_each_Vertex_within(this,vertex,r){

		buffer_fast_Read_Vector3D(recvBuffer,vertex->A);
		buffer_fast_Read_Vector3D(recvBuffer,vertex->Y);

	}end_for_each_Vertex_within

}

