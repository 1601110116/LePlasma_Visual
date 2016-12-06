/*
 * Visualize.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: zlstudio
 */

#include <Event.h>
#include <GL/freeglut_std.h>
#include <HostGrid.h>
#include <KeyBoard.h>
#include <omp.h>
#include <pthread.h>
#include <stddef.h>
#include <Vector3D.h>
#include <Visualize.h>
#include <cmath>
#include <sstream>
#include <vector>
#include <LePlasma.h>

//if cal speed is confined to 60 Hz
#define V_SYNC true

#define GridViewSizeX 500
#define GridViewSizeY 500

#define ExViewW 470
#define ExViewH 350

#define ExDistW 470
#define ExDistH 350

//distrubute
double E_scale=0.00005;
vector<double> ExDist;

int Exdx=0;
double Ex_size=0.00005;
double dt_scale=1;
vector<double> ExRecord;
vector<double> ExShow;
double ExMax=0;


pthread_t refresh_thread_handler;
pthread_mutex_t R_mutex,S_mutex;

bool _isHostView;
HostGrid* _Hostgrid;
Grid* _aimGrid;
void (*_idleFn)();
double fps;

double startTime;
int endTime;

bool running=true;

#include <CSPIC.h>


void* refresh_thread(void *p){
	while(true){
		pthread_mutex_lock(&R_mutex);
		double t=omp_get_wtime();
		if(running){
			if(_isHostView){
				_Hostgrid->askReport();

			}else if(_idleFn!=NULL){

				_idleFn();

			}else{
				pthread_mutex_unlock(&R_mutex);
				continue;
			}

			fps=floor(10/(omp_get_wtime()-t))/10;

			////Show

			double Ex=0;
			//line exists at x,5,5
			Ex=_aimGrid->Y_indexer[0][0][0]->x;

			ExRecord.push_back(Ex);

			//ExShow.push_back(log10(fabs(Ex/ExRecord[0])));
			ExShow.push_back(Ex);

			ExDist.clear();

			for (int i = 0; i < _aimGrid->gridX(); ++i) {
				ExDist.push_back(_aimGrid->Y_indexer[i][0][0]->x);
			}

		}
		/////
#if V_SYNC
		pthread_mutex_unlock(&S_mutex);
#else
		pthread_mutex_unlock(&R_mutex);
#endif
	}
}


Visualize::Visualize(Grid* grid,int ID):Stage("Node",320*ID+10,200,320,200),rx(0),ry(180),
		GridViewPort(),
		ExPlot(ExViewW,ExViewH,"x-asix","y-asix"),
		gridViewer(grid),
		coordinateIndicator(),
		aimGrid(grid),
		Hostgrid(NULL),isHostView(false),idleFn(NULL)
{

	gridViewer.x=320/2;
	gridViewer.y=200/2;
	gridViewer.viewTheta=180;
	//Temporary range set
	GridViewPort.width=300;
	GridViewPort.height=300;

	addChild(gridViewer);
	gridViewer.refreshGrid();

	pthread_mutex_init(&R_mutex,NULL);
	pthread_mutex_init(&S_mutex,NULL);
	pthread_mutex_lock(&S_mutex);

	_isHostView=isHostView;
	_Hostgrid=Hostgrid;
	_aimGrid=aimGrid;

	addEventListener(ENTER_FRAME,&EventListener::enterFrame_evt);
	addEventListener(MOUSE_WHEEL,&EventListener::mouseWheel_evt);
	addEventListener(MOUSE_DOWN,&EventListener::mouseDown_evt);
	addEventListener(MOUSE_DRAG,&EventListener::mouseDrag_evt);
	addEventListener(MOUSE_UP,&EventListener::mouseUp_evt);
	addEventListener(KEY_DOWN,&EventListener::keyDown_evt);


	pthread_create(&refresh_thread_handler,NULL,&refresh_thread,NULL);

}

Visualize::Visualize(Grid* grid):Stage("Node",400,300,1100,770),rx(0),ry(180),
		GridViewPort(),
		ExPlot(ExViewW,ExViewH,"tick","E"),
		ExDistribute(ExDistW,ExDistH,"x","y"),
		gridViewer(grid),
		coordinateIndicator(),
		aimGrid(grid),Hostgrid(NULL),
		idleFn(NULL),isHostView(false)
{
	initAnalysisView();
}

Visualize::Visualize(HostGrid* grid):Stage("Node",400,300,1000,770),rx(0),ry(180),
		GridViewPort(),
		ExPlot(ExViewW,ExViewH,"tick","E"),
		ExDistribute(ExDistW,ExDistH,"x","E"),
		gridViewer(grid),
		coordinateIndicator(),
		aimGrid(grid),Hostgrid(grid),
		idleFn(NULL),isHostView(true)
{

	initAnalysisView();
}

void Visualize::initAnalysisView(){
	GridViewPort.x=10;
	ExPlot.y=GridViewPort.y=25;
	ExDistribute.y=ExViewH+50;

	GridViewPort.graphics.lineStyle(0,0x0099ff,1);
	GridViewPort.graphics.drawRect(0,0,GridViewSizeX,GridViewSizeY);

	ExDistribute.x=ExPlot.x=40+GridViewSizeX;
	ExPlot.setCoordinate(PlotRange{Exdx*dt_scale,(Exdx+ExViewW)*dt_scale,-Ex_size,Ex_size});
	ExDistribute.setCoordinate(PlotRange{0,aimGrid->gridX(),-E_scale,E_scale});

	//ExDistribute.setCoordinate(PlotRange{-3,3,-E_scale,E_scale});
	//ExDistribute.drawFunction(_W1,0xff9900);
	//ExDistribute.drawFunction(dW1,0xff0099);


	gridViewer.x=GridViewSizeX/2+10;
	gridViewer.y=GridViewSizeY/2+12;

	coordinateIndicator.x=40;
	coordinateIndicator.y=GridViewSizeY-20;
	coordinateIndicator.viewTheta=gridViewer.viewTheta=180;

	//Temporary range set
	GridViewPort.width=GridViewSizeX;
	GridViewPort.height=GridViewSizeY;

	details.x=10;
	details.y=GridViewSizeY+50;
	details.setTextStyle(GLUT_BITMAP_HELVETICA_12,0x208800);

	addChild(gridViewer);
	addChild(coordinateIndicator);
	addChild(GridViewPort);
	addChild(ExPlot);
	addChild(ExDistribute);
	addChild(details);

	gridViewer.refreshGrid();

	pthread_mutex_init(&R_mutex,NULL);
	pthread_mutex_init(&S_mutex,NULL);
	pthread_mutex_lock(&S_mutex);

	_isHostView=isHostView;
	_Hostgrid=Hostgrid;
	_aimGrid=aimGrid;

	pthread_create(&refresh_thread_handler,NULL,&refresh_thread,NULL);


	addEventListener(ENTER_FRAME,&EventListener::enterFrame_evt);
	addEventListener(MOUSE_WHEEL,&EventListener::mouseWheel_evt);
	addEventListener(MOUSE_DOWN,&EventListener::mouseDown_evt);
	addEventListener(MOUSE_DRAG,&EventListener::mouseDrag_evt);
	addEventListener(MOUSE_UP,&EventListener::mouseUp_evt);
	addEventListener(KEY_DOWN,&EventListener::keyDown_evt);

}

Visualize::~Visualize() {
	// TODO Auto-generated destructor stub
}

void Visualize::enterFrame_evt(const Event&){

	_idleFn=idleFn;

#if V_SYNC
	int chk=pthread_mutex_trylock(&S_mutex);
#else
	int chk=pthread_mutex_trylock(&R_mutex);
#endif

	if(chk==0){

		gridViewer.refreshGrid();

		pthread_mutex_unlock(&R_mutex);
	}

	graphics.clear();
	graphics.textStyle(GLUT_BITMAP_HELVETICA_18,0xffffff);
	graphics.printAt(GridViewSizeX/2-35,5,"Grid Status");

	if(ExShow.size()>0){
		ExPlot.trackData1D(ExShow);

		endTime=omp_get_wtime()-startTime;
		int hours=endTime/3600;
		int mins=(endTime%3600)/60;
		int seconds=endTime%60;

		ostringstream output;
		output<<"run "<<ExShow.size()<<",  fps:"<<fps<<endl
				<<"Ex="<<ExRecord[ExRecord.size()-1]<<",  E0="<<ExRecord[0]<<endl
				<<"Time pass: "<<hours<<":"<<mins<<":"<<seconds;

		details.setText(output.str());
	}else{
		details.setText("Waiting for first frame...");
		startTime=omp_get_wtime();
	}
	if(ExDist.size()>0){
		ExDistribute.trackData1D(ExDist);
	}

}

void Visualize::resize_evt(const WindowEvent&){

}

bool GridDragging=false;
bool ExDragging=false;
double Ori_Ex_size;
int srx,sry;

void Visualize::mouseDrag_evt(const MouseEvent& evt){
	if(GridDragging){
		rx=srx+evt.x-evt.sx;
		ry=sry+evt.y-evt.sy;
		if(ry<=90)ry=90;
		if(ry>=270)ry=270;
		coordinateIndicator.viewPhi=gridViewer.viewPhi=rx;
		coordinateIndicator.viewTheta=gridViewer.viewTheta=ry;
	}else if(ExDragging){
		Exdx=srx+(evt.sx-evt.x);
		if(Exdx<0)Exdx=0;

		int dy=evt.sy-evt.y;
		Ex_size=Ori_Ex_size*pow(0.99,dy);

		ExPlot.setCoordinate(PlotRange{Exdx*dt_scale,(Exdx+ExViewW)*dt_scale,-Ex_size,Ex_size});
	}
}
void Visualize::mouseDown_evt(const MouseEvent& evt){
	if(GridViewPort.hitTest(evt.x,evt.y)){
		srx=rx;
		sry=ry;
		GridDragging=true;
	}else if(ExPlot.hitTest(evt.x,evt.y)){
		srx=Exdx;
		Ori_Ex_size=Ex_size;
		ExDragging=true;
	}
}

void Visualize::mouseUp_evt(const MouseEvent& evt){
	ExDragging=false;
	GridDragging=false;
}

void Visualize::mouseWheel_evt(const MouseEvent& evt){
	double delta=1+0.05*evt.delta;
	if(GridViewPort.hitTest(evt.x,evt.y)){
		gridViewer.cellSize*=delta;
		gridViewer.refreshGrid();
	}else if(ExPlot.hitTest(evt.x,evt.y)){
		dt_scale*=delta;
		ExPlot.setCoordinate(PlotRange{Exdx*dt_scale,(Exdx+ExViewW)*dt_scale,-Ex_size,Ex_size});
	}else if(ExDistribute.hitTest(evt.x,evt.y)){
		E_scale*=delta;
		ExDistribute.setCoordinate(PlotRange{0,aimGrid->gridX(),-E_scale,E_scale});

		//ExDistribute.setCoordinate(PlotRange{-3,3,-E_scale,E_scale});
	}
}
void Visualize::keyDown_evt(const KeyBoardEvent& evt){
	switch(evt.key){
		case LEFT:
			rx+=5;
			break;
		case RIGHT:
			rx-=5;
			break;
		case UP:
			ry-=5;
			break;
		case DOWN:
			ry+=5;
			break;
		case SPACE:
			running=!running;

	}
	rx=rx/5*5;
	ry=ry/5*5;
	if(ry<=90)ry=90;
	if(ry>=270)ry=270;
	coordinateIndicator.viewPhi=gridViewer.viewPhi=rx;
	coordinateIndicator.viewTheta=gridViewer.viewTheta=ry;
}
