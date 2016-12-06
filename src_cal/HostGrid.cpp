/*
 * HostGrid.cpp
 *
 *  Created on: Apr 13, 2016
 *      Author: zlstudio
 */

#include <Cell.h>
#include <HostGrid.h>
#if REPORT
#include <inet.h>
#endif
#include <Vertex.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


static ofstream outputFile;

HostGrid::HostGrid(int x,int y,int z,int n):Grid(x,y,z,true),nodes(n),workLength(x/n){

#if REPORT

	reportSocket=connectSocket(TRANSMITTER_ADDRESS,TRANSMITTER_OUT_PORT);
	if(reportSocket>0)	cout<<"Connected to Transmitter"<<endl;
#endif

	space=new Range[n]();

	for (int i = 0; i < n; ++i) {
		space[i]=Range(workLength*i,workLength*(i+1),0,_height,0,_length);
	}

#if RECORD
	time_t t = time(0);
	char tmp[64];
	strftime(tmp, sizeof(tmp), "%Y_%m_%d_%H_%M",localtime(&t));
	ostringstream name;
	name<<"Report "<<tmp<<".output";
	outputFile.open(name.str().c_str());
#endif

}

HostGrid::~HostGrid() {
}

void HostGrid::askReport(){
#if REPORT

	//	int cmd=9;
	//	send(reportSocket,&cmd,4,0);

	Cell* c;
	for_each_Cell(this,c){

		_particle_amount-=c->length();

		c->destroyAllParticlesInside();

	}end_for_each_Cell;

	for (int i = 0; i < nodes; ++i) {

		//get incoming Node ID
		int src;

		if(!Recv(reportSocket,&src,4))return;

		//cout<<"receiving.."<<endl;

		receiveParticlesFrom(src);
		receiveVertexesFrom(src);

	}

#endif

#if RECORD

	Vertex* v;

	//for_each_Vertex_within(grid,v,grid->workSpace){
	for_each_Vertex(this,v){

		outputFile<<v->Y.x<<'\t';

	}end_for_each_Vertex_within

	outputFile<<'\n';

#endif

}


bool HostGrid::receiveVertexesFrom(int src){

#if REPORT
	Vertex* v;

	int n;

	if(!Recv(reportSocket,&n,4))return false;

	if(n!=space[src].Volumn()){
		cout<<"vertexes do not match!"<<endl;
		return false;
	}

	//cout<<"read vertex buffer"<<endl;

	double recvbuffer[6];

	for_each_Vertex_within(this,v,space[src]){

		buffer_init();

		if(!Recv(reportSocket,recvbuffer,48))return false;//6*sizeof(double)

		buffer_Read_Vector3D(recvbuffer,v->A);
		buffer_Read_Vector3D(recvbuffer,v->Y);

	}end_for_each_Vertex

#endif
	return true;
}

bool HostGrid::receiveParticlesFrom(int src){

#if REPORT
	Electron* e_in;

	int n;

	if(!Recv(reportSocket,&n,4))return false;

	double recvbuffer[6];

	//cout<<"read particle buffer"<<endl;

	for (int j = 0; j < n; ++j) {

		buffer_init();

		if(!Recv(reportSocket,recvbuffer,48))return false;//6*sizeof(double)
		e_in=new Electron();

		buffer_Read_Vector3D(recvbuffer,e_in->Position);
		buffer_Read_Vector3D(recvbuffer,e_in->Momentum);

		//modify the position to global coordinate here.
		e_in->Position.x+=src*workLength-SWAP_LEN;

		//cout<<src<<" : add "<<e_in->Position.toString()<<endl;

		addParticle(e_in);
	}
	//cout<<"get"<<n<<endl;

#endif
	return true;
}
