/*
 * Event.h
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#ifndef INCLUDE_EVENT_H_
#define INCLUDE_EVENT_H_

typedef int eventType;

///MouseEvents:
#define MOUSE_DOWN 0
#define MOUSE_UP 1
#define MOUSE_DRAG 2

#define RIGHT_MOUSE_DOWN 3
#define RIGHT_MOUSE_UP 4
#define RIGHT_MOUSE_DRAG 5

#define MIDDLE_MOUSE_DOWN 6
#define MIDDLE_MOUSE_UP 7
#define MIDDLE_MOUSE_DRAG 8

#define MOUSE_MOVE 9

#define MOUSE_WHEEL 10


///KeyBoardEvents:
#define KEY_DOWN 11
#define KEY_UP 12


///Event:
#define ENTER_FRAME 13
#define RESIZE 14


#define EVENT_TYPES 15

struct MouseEvent{
		int x,y; //for point operations
		int sx,sy; //for drag operatons
		int delta; //for wheel
};

struct KeyBoardEvent{
		unsigned char key;//keyCode
		int x,y;//current mouseLocation

};

struct WindowEvent{
		int width,height;
};

struct Event{
		int tag;
};


#endif /* INCLUDE_EVENT_H_ */
