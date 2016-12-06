/*
 * Case.cpp
 *
 *  Created on: 2016-5-12
 *      Author: zlstudio
 */

#include "Grid.h"
#include <ctime>
#include "Case.h"
#include "RunManager.h"


Case::Case(){
	// TODO Auto-generated constructor stub

}

Case::~Case() {
	// TODO Auto-generated destructor stub
}

void Case::distributeParticle(){}
void Case::initP(){}
void Case::initA(){}
void Case::initY(){}

void Case::report(){}

void Case::launch(bool report){

	srand(0);

	if(RunManager::MPI_ID==0){
		distributeParticle();
	}

	grid->refreshParticleLocation();

	cout<<"Node "<<RunManager::MPI_ID<<" Distribute Complete. Particles: "<<grid->particles()<<endl;

	initA();
	initP();
	initY();

	if(RunManager::Nodes>1){
		grid->refreshParticleLocation();
	}

	if(report){
		time_t t = time(0);
		char tmp[64];
		strftime(tmp, sizeof(tmp), "%Y_%m_%d_%H_%M",localtime(&t));
		ostringstream name;
		name<<"Report "<<tmp<<".output";
		outputFile.open(name.str().c_str());
	}

}
