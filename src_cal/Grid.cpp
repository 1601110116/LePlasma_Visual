/*
 * Grid.cpp
 *
 *  Created on: Mar 3, 2016
 *      Author: zlstudio
 */

#include <Grid.h>

Grid::Grid(int x,int y,int z, bool period):workSpace(0,x,0,y,0,z),World(0,x,0,y,0,z){
	// TODO Auto-generated constructor stub
	initGrid3D(x,y,z);
	_period=period;

	particleSwap=new Cell();

	if(period)cout<<"Period ";
	cout<<"Grid initialized, size: "<<x<<"*"<<y<<"*"<<z<<endl;

	//TODO: temp scale

	_workLength=x;
	_scale=2.4e-4*METER;
}

Grid::~Grid() {
}

void Grid::initGrid3D(int x,int y,int z){

	VertexContainer=new Vertex***[x];
	A_indexer=new Vector3D***[x];
	Y_indexer=new Vector3D***[x];

	for (int i = 0; i < x; ++i) {
		Vertex*** d2=new Vertex**[y];
		Vector3D*** d2_A=new Vector3D**[y];
		Vector3D*** d2_Y=new Vector3D**[y];
		for (int j = 0; j < y; ++j) {
			Vertex** d1=new Vertex*[z];
			Vector3D** d1_A=new Vector3D*[z];
			Vector3D** d1_Y=new Vector3D*[z];
			for (int k = 0; k < z; ++k) {
				Vertex* d0=new Vertex(i,j,k);
				d1[k]=d0;
				d1_A[k]=&d0->A;
				d1_Y[k]=&d0->Y;
			}
			d2[j]=d1;
			d2_A[j]=d1_A;
			d2_Y[j]=d1_Y;
		}
		VertexContainer[i]=d2;
		A_indexer[i]=d2_A;
		Y_indexer[i]=d2_Y;
	}

	//priod condition

	CellContainer=new Cell***[x];

	for (int i = 0; i < x; ++i) {
		Cell*** d2=new Cell**[y];
		for (int j = 0; j < y; ++j) {
			Cell** d1=new Cell*[z];
			for (int k = 0; k < z; ++k) {
				Cell* d0=new Cell(i,j,k);
				d1[k]=d0;
			}
			d2[j]=d1;
		}
		CellContainer[i]=d2;

	}

	_width=x;
	_height=y;
	_length=z;
	_particle_amount=0;
}

bool Grid::addParticle(Particle* p){


	int cell_x=p->Position.x;
	int cell_y=p->Position.y;
	int cell_z=p->Position.z;

	if(CellContainer[cell_x][cell_y][cell_z]->pushParticle(p)){
		_particle_amount++;
		return true;
	}


	//No longer test position

	//	int cell_x=floor(p->Position.x);
	//	int cell_y=floor(p->Position.y);
	//	int cell_z=floor(p->Position.z);



	//	if(_period){
	//		p->Position.x-=getPeriod(cell_x,_width)*_width;
	//		p->Position.y-=getPeriod(cell_y,_height)*_height;
	//		p->Position.z-=getPeriod(cell_z,_length)*_length;
	//	}

	//	if(cell(cell_x,cell_y,cell_z)->pushParticle(p)){
	//		_particle_amount++;
	//
	//		return true;
	//	}

	cerr<<"Warning from Grid : Adding Particle Failed!"<<endl;

	return false;
}

bool Grid::directAddParticle(Particle* p){

	if(CellContainer[0][0][0]->pushParticle(p)){
		_particle_amount++;
		return true;
	}

	cerr<<"Warning from Grid : Adding Particle Failed!"<<endl;

	return false;
}

bool Grid::moveParticleTo(Particle* p, const Vector3D &location){
	if(p->cell!=NULL){

		int cell_x=floor(p->Position.x);
		int cell_y=floor(p->Position.y);
		int cell_z=floor(p->Position.z);

		p->Position=location;

		if(_period){
			p->Position.x-=getPeriod(cell_x,_width)*_width;
			p->Position.y-=getPeriod(cell_y,_height)*_height;
			p->Position.z-=getPeriod(cell_z,_length)*_length;
		}

		if(cell(cell_x,cell_y,cell_z)==p->cell){
			return true;
		}


		if(!p->cell->withdrawParticle(p)){
			cerr<<"Warning from Grid : Moving Particle Failed!"<<endl;
			return false;
		}

		cell(cell_x,cell_y,cell_z)->pushParticle(p);

		return true;
	}
	cerr<<"Warning from Grid : Moving Particle Failed!"<<endl;

	return false;
}

void Grid::refreshParticleLocation(){


	Particle* p,*p_swap;
	for_each_Particle(this,p){

		if(p->cell!=NULL){

			int cell_x=floor(p->Position.x);
			int cell_y=floor(p->Position.y);
			int cell_z=floor(p->Position.z);

			if(_period){
				p->Position.x-=getPeriod(cell_x,_width)*_width;
				p->Position.y-=getPeriod(cell_y,_height)*_height;
				p->Position.z-=getPeriod(cell_z,_length)*_length;
			}

			if(cell(cell_x,cell_y,cell_z)!=p->cell){

				p_swap=p;

				p=p->nextParticle;

				//we have to push this particle into Swap first in order to complete the iteration, or some of the particles will be ignored.
				if(!p_swap->cell->withdrawParticle(p_swap)){
					cerr<<"Warning from Grid : Withdrawing Particle Failed!"<<endl;
				}

				if(!particleSwap->pushParticle(p_swap)){
					cerr<<"Warning from Grid : Pushing Particle to Swap Failed!"<<endl;
				}
			}else{
				p=p->nextParticle;
			}
		}else{
			cerr<<"Warning from Grid : Particle with no cell pointer exists!"<<endl;
		}

	}end_for_each_Particle_raw

	//put the swap particles to where they should go.

	for_each_Particle_in_Cell(particleSwap,p){

		p_swap=p;

		p=p->nextParticle;

		particleSwap->withdrawParticle(p_swap);
		cell(p_swap->Position.x,p_swap->Position.y,p_swap->Position.z)->pushParticle(p_swap);

	}end_for_each_Particle_in_Cell_raw;

}

void Grid::syncSharedParticles(){}
void Grid::reportToHost(){}

int Grid::particles(){
	return _particle_amount;
}

//TODO: %must die

Cell* Grid::cell(int a,int b ,int c){

#if OPTIMIZE_1D_MPI
	b=c=0;
#elif OPTIMIZE_1D
	//As we can confine the rang in +/-2
	//for height=lenght=1 and width>=2
	a<0?(a+=_width):(a>=_width?a-=_width:a);
	b=c=0;
#elif OPTIMIZE_3D
	a<0?(a+=_width):(a>=_width?a-=_width:a);
	b<0?(b+=_height):(b>=_height?b-=_height:b);
	c<0?(c+=_length):(c>=_length?c-=_length:c);
#else
	if(_period){
		a=a%_width;
		b=b%_height;
		c=c%_length;

		if(a<0)a+=_width;
		if(b<0)b+=_height;
		if(c<0)c+=_length;
	}
#endif

#if DEGUB
	if(a<_width&&b<_height&&c<_length&&a>=0&&b>=0&&c>=0){
		return CellContainer[a][b][c];
	}else{
		cerr<<"Warning from Grid : Requesting unexisted Cell at ("<<a<<","<<b<<","<<c<<")"<<endl;
		return NULL;
	}
#else
	return CellContainer[a][b][c];
#endif

}
Vertex* Grid::vertex(int a, int b, int c){

#if OPTIMIZE_1D_MPI
	b=c=0;
#elif OPTIMIZE_1D
	//As we can confine the rang in +/-2
	//for height=lenght=1 and width>=2
	a<0?(a+=_width):(a>=_width?a-=_width:a);
	b=c=0;
#elif OPTIMIZE_3D
	a<0?(a+=_width):(a>=_width?a-=_width:a);
	b<0?(b+=_height):(b>=_height?b-=_height:b);
	c<0?(c+=_length):(c>=_length?c-=_length:c);
#else
	if(_period){
		a=a%_width;
		b=b%_height;
		c=c%_length;

		if(a<0)a+=_width;
		if(b<0)b+=_height;
		if(c<0)c+=_length;
	}
#endif


#if DEGUB
	if(a<_width&&b<_height&&c<_length&&a>=0&&b>=0&&c>=0){
		return VertexContainer[a][b][c];
	}else{
		cerr<<"Warning from Grid : Requesting unexisted Vertex at ("<<a<<","<<b<<","<<c<<")"<<endl;
		return NULL;
	}
#else
	return VertexContainer[a][b][c];
#endif

}

int Grid::gridX(){
	return _width;
}

int Grid::gridY(){
	return _height;
}
int Grid::gridZ(){
	return _length;
}

int Grid::lengthX(){
	return _width;
}
int Grid::lengthY(){
	return _height;
}
int Grid::lengthZ(){
	return _length;
}

double Grid::scale(){
	return _scale;
}

int Grid::workLength(){
	return _workLength;
}

void Grid::showGridMap(){
	cout<<endl<<"Grid Map: "<<_width<<"*"<<_height<<"*"<<_length<<endl<<"cells:"<<endl;
	for (int i = 0; i < _width; ++i) {
		for (int j = 0; j < _height; ++j) {
			for (int k = 0; k < _length; ++k) {
				Cell* _cell=cell(i,j,k);
				cout<<_cell->length()<<"\t";
			}
			cout<<endl;
		}
		cout<<endl;
	}
}
