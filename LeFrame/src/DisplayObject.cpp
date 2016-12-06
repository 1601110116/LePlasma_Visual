/*
 * DisplayObject.cpp
 *
 *  Created on: Apr 24, 2016
 *      Author: zlstudio
 */

#include <DisplayObject.h>
#include <Graphics.h>
#include <stddef.h>
#include <Stage.h>
#include <iterator>
#include <vector>

DisplayObject::DisplayObject():Rect(0,0,0,0),viewDistance(1),viewTheta(0),viewPhi(0){
	// TODO Auto-generated constructor stub

}

DisplayObject::~DisplayObject() {
	// TODO Auto-generated destructor stub
}

void DisplayObject::addChild(DisplayObject* child){
	childList.push_back(child);
}

void DisplayObject::addChild(DisplayObject& child){
	childList.push_back(&child);
}

void DisplayObject::removeChild(DisplayObject* child){
	child->stage=NULL;
	for (unsigned int i = 0; i < childList.size(); ++i) {
		if(childList[i]==child){
			childList.erase(childList.begin()+i);
		}
	}
}

void DisplayObject::removeChildAt(int i){
	childList[i]->stage=NULL;
	childList.erase(childList.begin()+i);

}

//cx: coordinate X
void DisplayObject::redraw(double cx,double cy){
	//first draw itself

	globalX=cx+x;
	globalY=cy+y;

	moveCoordinate(globalX,stage->stageHeight-globalY,0);
	seeFrom(viewDistance,viewTheta,viewPhi);

	graphics.draw();

	//then draw its children
	for (unsigned int i = childList.size(); i > 0; --i) {
		childList[i-1]->stage=stage;
		childList[i-1]->redraw(globalX,globalY);
	}
}

bool DisplayObject::hitTest(double _x,double _y){
	return (_x>globalX)&&_x<(width+globalX)&&_y>globalY&&_y<globalY+height;
}
