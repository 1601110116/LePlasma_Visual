/*
 * Buffer.h
 *
 *  Created on: Apr 14, 2016
 *      Author: zlstudio
 */
#include "mpi.h"
#ifndef INCLUDE_BUFFER_H_
#define INCLUDE_BUFFER_H_


#define length_of(array) sizeof(array)/sizeof(array[0])

#define buffer_fast_init(name,length)\
		int buf_n_##name=-1;\
		double name[length];

#define buffer_init()\
		int buffer_n=-1;\

#define buffer_Write_Vector3D(buffer,vector)\
		buffer[++buffer_n]=vector.x;\
		buffer[++buffer_n]=vector.y;\
		buffer[++buffer_n]=vector.z;\

#define buffer_Read_Vector3D(buffer,vector)\
		vector.x=buffer[++buffer_n];\
		vector.y=buffer[++buffer_n];\
		vector.z=buffer[++buffer_n];\


#define buffer_fast_Write_Vector3D(buffer,vector)\
		buffer[++buf_n_##buffer]=vector.x;\
		buffer[++buf_n_##buffer]=vector.y;\
		buffer[++buf_n_##buffer]=vector.z;\

#define buffer_fast_Read_Vector3D(buffer,vector)\
		vector.x=buffer[++buf_n_##buffer];\
		vector.y=buffer[++buf_n_##buffer];\
		vector.z=buffer[++buf_n_##buffer];\

#define buffer_fast_Write_Double(buffer,a)\
		buffer[++buf_n_##buffer]=a;\

#define buffer_fast_Read_Double(buffer) buffer[++buf_n_##buffer];


#define buffer_fast_send_to(buffer, dest,tag) MPI_Send(buffer,length_of(buffer),MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
#define buffer_send_to(buffer,len,dest,tag) MPI_Send(buffer,len,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);

#define buffer_fast_recv_from(buffer, src, tag)\
		MPI_Status Dummystate;\
		MPI_Recv(buffer,length_of(buffer),MPI_DOUBLE,src,tag,MPI_COMM_WORLD,&Dummystate);


#endif /* INCLUDE_BUFFER_H_ */
