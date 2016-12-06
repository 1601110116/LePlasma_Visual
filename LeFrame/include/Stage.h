/*
 * Stage.h
 *
 *  Created on: Apr 24, 2016
 *      Author: zlstudio
 */

#ifndef LEFRAME_INCLUDE_STAGE_H_
#define LEFRAME_INCLUDE_STAGE_H_

#include <Sprite.h>

class Stage:public Sprite {
	public:
		Stage();

		Stage(const char* title,int x,int y,int w,int h);

		virtual ~Stage();

		void display();

		const char* title;

		Stage *root;

		int stageX,stageY;
		int stageHeight,stageWidth;

};

#endif /* LEFRAME_INCLUDE_STAGE_H_ */
