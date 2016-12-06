/*
 * GridManager.h
 *
 *  Created on: Apr 6, 2016
 *      Author: zlstudio
 */

#ifndef INCLUDE_MPIGRID_H_
#define INCLUDE_MPIGRID_H_



#if REPORT
#include <inet.h>
#endif



#include <mpi.h>
#include "Grid.h"

#define MPItrace(inf) cout<<"[Node "<<thisNode<<"] "<<inf<<endl;fflush(stdout);

#define SYNC_VERTEX 100

#define SYNC_PARTICLE 200

#define SHARE_PARTICLE 250

#define SYNC_CACHE  300

#define SIGNAL 300

#define SWAP_LEN 3

#define V_ADJ_BEGIN 0
#define V_SHARE_END 6

#define P_ADJ_BEGIN_L 1
#define P_ADJ_BEGIN_R 2
#define P_SHARE_END_L 4
#define P_SHARE_END_R 5

////////MPI_GRID architecture
////////prevNode <-- thisNode --> nextNode
////////Architecture for each Node:

////////|0	|1   	|2 	|3 	|4	...	...	... 	 ...	...		|-5	|-4	|-3	|-2	|-1 	|///////worklength+6
////////|	   	 	 	 			storeLength 	 			 		 	|///////worklength+6
////////|      SWAP	|				workSpace				|     SWAP	|///////workSpace&SWAP
////////|       VadjL	|        VshL	|				|       VshR	|       VadjR	|///////Vertex shared(sh)/adjacent(adj) space
////////|	|   PadjL	|PshL|							|   PshR	|PadjR|		|///////Particle shared(sh)/adjacent(adj) space
////////| SWAP_LEN	|<==SWAP_LEN	storeLength-SWAP_LEN==>	| SWAP_LEN	|


inline int wait_signal_from(int src){
	int signal;
	MPI_Status Dummystate;
	MPI_Recv(&signal,1,MPI_INT,src,SIGNAL,MPI_COMM_WORLD,&Dummystate);
	return signal;
}

inline void send_signal_to(int dest, int signal){
	MPI_Send(&signal,1,MPI_INT,dest,SIGNAL,MPI_COMM_WORLD);
}


using namespace std;


class Distributer;

class MPIGrid : public Grid{
	public:
		MPIGrid(int w,int h,int l);
		virtual ~MPIGrid();


		//override:
		bool directAddParticle(Particle* p);
		void refreshParticleLocation();
		void syncSharedParticles();


		//for Report, Now we consider to possibilities. if socket connection is possible, we can run a host on laptop. Or we will report to Node 0.
		void reportToHost();
		void report();

	protected:

		//for A & Y
		void syncVertex();

		//for particles (X & P)
		void syncParticles();

		//for share
		void shareParticles();

		Distributer* distributer;

		///send vertexes within range to certain destination
		void sendVertexes(int dest,const Range&);

		///receive vertexes from certain source to refresh certain range
		void recvVertexes(int source, const Range&);

		///send particles within range to certain destination
		void sendParticles(int dest,const Range&);

		///receive particles from certain source to refresh certain range
		void recvParticles(int source,const Range&);


		//for particles Momentum
		void sendMomentum(int dest,const Range&);
		void recvMomentum(int source,const Range&);

		///for report
		void reportParticles();
		void reportVertexes();

#if REPORT
		int report_socket;
#endif

		int nodes;
		int thisNode;
		int nextNode;
		int prevNode;

		int _worldwidth,_worldheight,_worldlength;

		int _storeLength;

		Range VxAdjacentL,VxSharedL,VxSharedR,VxAdjacentR;
		Range PtAdjacentL,PtSharedL,PtSharedR,PtAdjacentR;

		friend class ParticleDistributer;
		friend class CaseLauncher;
		friend class MPIRunManager;
};

#endif /* INCLUDE_MPIGRID_H_ */
