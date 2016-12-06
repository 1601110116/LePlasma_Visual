/*
 * Plot.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: zlstudio
 */

#include <Plot.h>

Plot::Plot() {
	// TODO Auto-generated constructor stub
}

Plot::~Plot() {
	// TODO Auto-generated destructor stub
}

Plot::Plot(int w,int h,const char* xlab ,const char* ylab):baseX(0),baseY(0),basescaleX(1),basescaleY(1)
{

	xname=xlab;
	yname=ylab;

	//		addEventListener(Event.ADDED_TO_STAGE,function (e):void{
	//			stage.addEventListener(MouseEvent.MOUSE_MOVE,function (e):void{
	//				mouseTracker.setXY(Mathx::cutFloat((mouseX/xscale+SX)/basescaleX,4),Mathx::cutFloat(((Height-mouseY)/yscale+SY)/basescaleY,4));
	//			});
	//		});

	autoRange=false;

	addChild(graph);

	setSize(w,h);
}

void Plot::setCoordinate(const PlotRange &range){

	if(range.sx==range.ex||range.sy==range.ey)return;

	graphics.clear();
	graphics.lineStyle(0,FRAME_COLOR,0.2);
	graphics.drawRect(0,0,Width,Height);

	double tickX,tickY;
	double tickSX,tickEX;
	double tickSY,tickEY;

	if(autoRange){

		Arrange arrangeX=Mathx::Tick125(range.sx,range.ex);
		tickSX=SX=arrangeX.start;
		tickEX=EX=arrangeX.end;
		tickX=arrangeX.tick;

		Arrange arrangeY=Mathx::Tick125(range.sy,range.ey);
		tickSY=SY=arrangeY.start;
		tickEY=EY=arrangeY.end;
		tickY=arrangeY.tick;

	}else{
		SX=range.sx;
		EX=range.ex;
		SY=range.sy;
		EY=range.ey;
		Arrange arrangeX=Mathx::Tick125(range.sx,range.ex);
		Arrange arrangeY=Mathx::Tick125(range.sy,range.ey);

		tickX=arrangeX.tick;
		tickY=arrangeY.tick;

		tickSX=arrangeX.start+tickX;
		tickEX=arrangeX.end-tickX;
		tickSY=arrangeY.start+tickY;
		tickEY=arrangeY.end-tickY;

		if(tickX==0||tickY==0){
			cerr<<"Tick equals 0!!"<<endl;return;
		}

	}

	baseX=floor(log10(range.ex-range.sx));
	if(baseX>=4||baseX<=-4){
		xlabel=xname+"*10^"+Mathx::DoubleToString(baseX);
		basescaleX=pow(10,baseX);
	}else{
		xlabel=xname;
		basescaleX=1;
	}

	baseY=floor(log10(range.ey-range.sy));
	if(baseY>=4||baseY<=-4){
		ylabel=yname+"*10^"+Mathx::DoubleToString(baseY);
		basescaleY=pow(10,baseY);
	}else{
		ylabel=yname;
		basescaleY=1;
	}


	xscale=Width/(EX-SX);
	yscale=Height/(EY-SY);

	///Draw now
	graphics.textStyle(FONT_COLOR);
	graphics.lineStyle(0,AXIS_COLOR,0.6);

	//y-axis
	double axisX=Max(-SX*xscale,0);
	axisX=Min(axisX,Width);

	graphics.moveTo(axisX,0);
	graphics.lineTo(axisX,Height);
	//tick
	for (double i = tickSY; i <=tickEY ; i+=tickY){
		graphics.moveTo(axisX,Height-(i-SY)*yscale);
		graphics.lineTo(axisX+4,Height-(i-SY)*yscale);
		if(fabs(i)>=1e-10){
			graphics.printAt(axisX+4,Height-(i-SY)*yscale-8,Mathx::showFloat(i/basescaleY,3));
		}
	}

	//x-axis
	double axisY=Max(Height+SY*yscale,0);
	axisY=Min(axisY,Height);

	graphics.moveTo(0,axisY);
	graphics.lineTo(Width,axisY);

	for (double i = tickSX; i <tickEX ;i+=tickX){
		graphics.moveTo((i-SX)*xscale,axisY);
		graphics.lineTo((i-SX)*xscale,axisY-4);

		graphics.printAt((i-SX)*xscale,axisY,Mathx::showFloat(i/basescaleX,3));
	}

	///Axis-name
	graphics.printAt(0,-14,ylabel);
	graphics.printAt(Width/2,Height+5,xlabel);
	//redraw();
}
void Plot::setSize(int w,int h){
	Width=w;
	Height=h;

	xscale=Width/(EX-SX);
	yscale=Height/(EY-SY);

	width=Width;
	height=Height;
	setCoordinate(range);
}

void Plot::trackData1D(const vector<double>& data){

	graph.graphics.clear();

	unsigned int sx=Max(0,SX);
	unsigned int ex=Min(data.size(),EX);

	if(sx>data.size())return;

	graph.graphics.lineStyle(0,GRAPH_COLOR,1);
	graph.graphics.moveTo(0,Height-yscale*(data[sx]-SY));

	for (unsigned int i = sx; i < ex; i++) {
		graph.graphics.lineTo((i-SX)*xscale,Height-yscale*(data[i]-SY));
	}
}
//
//public function addFuntion(fn:Function,name:String,lstyle:LineStyle=null):void{
//	if(lstyle==null){
//		lstyle=LineStyle.defaultLineStyle;
//	}
//	//Analyse Range
//	if(autoRange){
//		if(noRange){
//			rearrange=true;
//			noRange=false;
//			Range=[SX,EX,Number.MAX_VALUE,Number.MIN_VALUE];
//		}
//		var steplength:Number=(EX-SX)/Width;
//		var steps:int=(EX-SX)/steplength;
//		var fnres:Number;
//		var rearrange:Boolean=false;
//		for (var i:Number = 0; i < steps+1; i++) {
//			fnres=fn(i*steplength+SX);
//			if(!isNaN(fnres)){
//				if(fnres>Range[3]){Range[3]=fnres;rearrange=true}
//				if(fnres<Range[2]){Range[2]=fnres;rearrange=true}
//			}
//		}
//		if(rearrange){
//			setCoordinate(Range,xtick,ytick);
//		}
//	}
//
//	legend.push(new FunctionGraph(fn,name,lstyle));
//	Items.push({type:ITEM_FUNCTION,fn:fn,name:name,lstyle:lstyle})
//													drawFunction(fn,name,lstyle);
//}
void Plot::drawFunction(double (*fn)(double x),Color c){

	double steplength=(EX-SX)/Width;
	int steps=(EX-SX)/steplength;

	double fnres;

	bool prevIsNaN=true;

	graphics.lineStyle(0,c,1);

	for (int i = 0; i < steps+1; i++) {
		fnres=fn(i*steplength+SX);
		if(isnan(fnres)){
			prevIsNaN=true;
		}else{
			if(prevIsNaN){
				graphics.moveTo(i,Height-yscale*(fnres-SY));
				prevIsNaN=false;
			}else{
				graphics.lineTo(i,Height-yscale*(fnres-SY));
			}
		}
	}
}


//}
