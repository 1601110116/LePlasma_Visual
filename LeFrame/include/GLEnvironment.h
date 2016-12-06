/*
 * GLEnvironment.h
 *
 *  Created on: Apr 24, 2016
 *      Author: zlstudio
 */

#ifndef OPENGL_GLENVIRONMENT_H_
#define OPENGL_GLENVIRONMENT_H_

#include <GL/glut.h>
#include <iostream>
#include <Stage.h>
#include <Event.h>
#include <GLEnvironment.h>
#include <GLEventManager.h>

using namespace std;

static Stage* _stage;

static void MouseHandler(int button,int state,int x,int y);
static void MouseDragHandler(int x,int y);
static void MouseMoveHandler(int x,int y);

static void KeyDownHandler(unsigned char key, int x, int y);
static void KeyUpHandler(unsigned char key, int x, int y);

static void SpecialKeyDownHandler(int key, int x, int y);
static void SpecialKeyUpHandler(int key, int x, int y);

static bool L_draging;
static bool R_draging;
static bool M_draging;

static int L_startDragX,L_startDragY;
static int R_startDragX,R_startDragY;
static int M_startDragX,M_startDragY;

void InitGLEnvironment();
void InitStage(Stage* stage);
static inline void _reshape(int w, int h);
static inline void _display();
static inline void _idle();

#endif /* OPENGL_GLENVIRONMENT_H_ */
