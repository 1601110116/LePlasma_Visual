/*
 * Sprite.h
 *
 *  Created on: Apr 23, 2016
 *      Author: zlstudio
 */

#ifndef OPENGL_SPRITE_H_
#define OPENGL_SPRITE_H_

#include <vector>
#include <iostream>
#include <DisplayObject.h>
#include <EventListener.h>

using namespace std;
class Sprite:public DisplayObject,public EventListener {
	public:
		Sprite();

		virtual ~Sprite();

	private:


};

#endif /* OPENGL_SPRITE_H_ */
