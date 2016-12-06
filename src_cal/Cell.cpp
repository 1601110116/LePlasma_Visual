/*
 * Cell.cpp
 *
 *  Created on: Mar 3, 2016
 *      Author: zlstudio
 *      miao~
 */

#include "Cell.h"

Cell::Cell():_x(0),_y(0),_z(0),head(NULL),_length(0),restorePointer(NULL),restoreLength(0),tail(NULL){
}

Cell::Cell(int x,int y,int z):_x(x),_y(y),_z(z),head(NULL),_length(0),restorePointer(NULL),restoreLength(0),tail(NULL) {
}

Cell::~Cell() {
	// TODO Auto-generated destructor stub
}

int Cell::x(){
	return _x;
}

int Cell::y(){
	return _y;
}

int Cell::z(){
	return _z;
}

string Cell::toString(){
	ostringstream tmpchar;

	tmpchar<<"Cell at ("<<_x<<","<<_y<<","<<_z<<")";

	return tmpchar.str();
}
Particle* Cell::firstParticle(){

	if(_length){
		return head;
	}
	return NULL;
}
#define SQUENCE true
bool Cell::pushParticle(Particle* p){ //insert to head  -- no longer insert to head as we now need to sync between nodes. Now we insert to tail (2016/5/11)

#if SQUENCE
	if(p->nextParticle==NULL&&p->prevParticle==NULL&&p->cell==NULL){

		if(_length>0){
			p->prevParticle=tail;
			tail->nextParticle=p;
		}else{
			head=p;
		}

		tail=p;


		p->cell=this;
		_length++;
		return true;
	}else{
		cerr<<"Warning from ParticleList : Attempting to push interdependent particle!"<<endl;
	}
	return false;
#else

	if(p->nextParticle==NULL&&p->prevParticle==NULL&&p->cell==NULL){

		if(_length>0){
			p->nextParticle=head;
			head->prevParticle=p;
		}
		head=p;

		p->cell=this;
		_length++;
		return true;
	}else{
		cerr<<"Warning from ParticleList : Attempting to push interdependent particle!"<<endl;
	}
	return false;

#endif
}

int Cell::overwriteParticle(Particle* p){

	if(restorePointer==NULL){
		//This means no more particle can be used
		Particle* tobeAdd=p->clone();
		tobeAdd->Position=p->Position;
		tobeAdd->Momentum=p->Momentum;
		pushParticle(tobeAdd);

		return 1;
	}else{
		restorePointer->Position=p->Position;
		restorePointer->Momentum=p->Momentum;
		restoreLength++;
		restorePointer=restorePointer->nextParticle;

		return 0;
	}
}

void Cell::restore(){
	restorePointer=head;
	restoreLength=0;
}
int Cell::cut(){
	if(restorePointer==NULL){
		//This means we do not need to cut tail(particles did not be overwrited).
		return 0;
	}else{

		//cut tail.
		if(restorePointer==head){
			head=NULL;
		}else{
			restorePointer->prevParticle->nextParticle=NULL;
		}

		Particle* iterator,*tobeDeleted;
		iterator=restorePointer;
		while(iterator!=NULL){
			tobeDeleted=iterator;
			iterator=iterator->nextParticle;

			delete tobeDeleted;

		}

		int dl=_length-restoreLength;
		_length=restoreLength;
		return dl;
	}
}

bool Cell::withdrawParticle(Particle* p){

#if SQUENCE
	if(p->cell==this){
		if(p==head){
			head=p->nextParticle;
			if(head!=NULL){
				head->prevParticle=NULL;
			}else{
				tail=NULL;
			}
		}else if(p==tail){
			p->prevParticle->nextParticle=NULL;
			tail=p->prevParticle;
		}else{
			p->prevParticle->nextParticle=p->nextParticle;
			p->nextParticle->prevParticle=p->prevParticle;
		}

		p->cell=NULL;
		p->prevParticle=NULL;
		p->nextParticle=NULL;

		_length--;
		return true;
	}else{
		cerr<<"Warning from ParticleList : Attempting to withdraw unexisted particle!"<<endl;
	}
	return false;
#else
	if(p->cell==this){
		if(p==head){
			head=p->nextParticle;
			if(head!=NULL){
				head->prevParticle=NULL;
			}
		}else if(p->nextParticle==NULL){
			p->prevParticle->nextParticle=NULL;
		}else{
			p->prevParticle->nextParticle=p->nextParticle;
			p->nextParticle->prevParticle=p->prevParticle;
		}

		p->cell=NULL;
		p->prevParticle=NULL;
		p->nextParticle=NULL;

		_length--;
		return true;
	}else{
		cerr<<"Warning from ParticleList : Attempting to withdraw unexisted particle!"<<endl;
	}
	return false;
#endif
}

void Cell::destroyAllParticlesInside(){
	Particle* iterator,*tobeDeleted;
	iterator=head;
	while(iterator!=NULL){
		tobeDeleted=iterator;
		iterator=iterator->nextParticle;

		delete tobeDeleted;

	}
	head=NULL;
	tail=NULL;
	_length=0;

}

int Cell::length(){
	return _length;
}


string Cell::treeView(){
	Particle* tmp=head;
	ostringstream tmpchar;

	tmpchar<<"Head:";

	while(tmp!=NULL){
		tmpchar<<tmp->name<<"->";
		tmp=tmp->nextParticle;
	}
	tmpchar<<"END"<<endl<<"Length:"<<_length<<endl;
	return  tmpchar.str();
}
