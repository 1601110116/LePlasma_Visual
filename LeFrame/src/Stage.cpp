/*
 * Stage.cpp
 *
 *  Created on: Apr 24, 2016
 *      Author: zlstudio
 */

#include <Stage.h>

Stage::Stage(const char* title,int x,int y,int w,int h):
title(title),
stageX(x),stageY(y),
stageWidth(w),stageHeight(h)
{
	this->x=0;
	this->y=0;
	stage=this;
}

Stage::Stage():title(title),
		stageX(0),stageY(0),
		stageWidth(DEFAULT_STAGE_WIDTH),
		stageHeight(DEFAULT_STAGE_HEIGHT)
{
	stage=this;
}

Stage::~Stage(){}

void Stage::display(){
	redraw(0,0);
}
