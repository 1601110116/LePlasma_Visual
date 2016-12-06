/*
 * GLEnvironment.cpp
 *
 *  Created on: Apr 24, 2016
 *      Author: zlstudio
 */

#include <GLEnvironment.h>

void InitGLEnvironment(){
	//////Central Event Generator & Dispatcher
	new GLEventManager();
}

void InitStage(Stage* stage){

	int argc=0;
	char* argv[]={};

	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);

	glutInitWindowPosition(stage->stageX, stage->stageY);

	glutInitWindowSize(stage->stageWidth, stage->stageHeight);

	glutCreateWindow(stage->title);


	///Environment Variables

	L_draging=false;
	R_draging=false;
	M_draging=false;

	L_startDragX=L_startDragY=R_startDragX=R_startDragY=M_startDragX=M_startDragY=0;
	_stage=stage;


	//Call backs
	glutReshapeFunc(_reshape);

//	cout<<"here"<<endl;

	glutDisplayFunc(_display);

	glutIdleFunc(_idle);


	glutMouseFunc(MouseHandler);
	glutMotionFunc(MouseDragHandler);
	glutPassiveMotionFunc(MouseMoveHandler);

	glutKeyboardFunc(KeyDownHandler);
	glutKeyboardUpFunc(KeyUpHandler);

	glutSpecialFunc(SpecialKeyDownHandler);
	glutSpecialUpFunc(SpecialKeyUpHandler);

	glutMainLoop();
}

void _reshape(int w, int h){
	_stage->stageWidth=w;
	_stage->stageHeight=h;

	glViewport(0, 0, (GLsizei)w, (GLsizei)h);

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

	if (w <= h)

		glOrtho(0, w, 0, h, -1000.0, 1000.0);

	else

		glOrtho(0, w, 0, h, -1000.0, 1000.0);

	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();

	GLEventManager::dispatcherList[RESIZE]->dispatch(WindowEvent{w,h});
}

void _display(){
}
void _idle(){

	GLEventManager::dispatcherList[ENTER_FRAME]->dispatch(Event{});

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	_stage->display();

	glFlush();

	glutSwapBuffers();
}


///////////////Event

void MouseHandler(int button,int state,int x,int y)
{

	int mouseState=button*10+state;

	switch(mouseState){
		case 0:
			L_startDragX=x;
			L_startDragY=y;
			L_draging=true;
			GLEventManager::dispatcherList[MOUSE_DOWN]->dispatch(MouseEvent{x,y});
			break;
		case 1:
			L_draging=false;
			GLEventManager::dispatcherList[MOUSE_UP]->dispatch(MouseEvent{x,y});
			break;
		case 10:
			M_startDragX=x;
			M_startDragY=y;
			M_draging=true;
			GLEventManager::dispatcherList[MIDDLE_MOUSE_DOWN]->dispatch(MouseEvent{x,y});
			break;
		case 11:
			M_draging=false;
			GLEventManager::dispatcherList[MIDDLE_MOUSE_UP]->dispatch(MouseEvent{x,y});
			break;
		case 20:
			R_startDragX=x;
			R_startDragY=y;
			R_draging=true;
			GLEventManager::dispatcherList[RIGHT_MOUSE_DOWN]->dispatch(MouseEvent{x,y});

			break;
		case 21:
			R_draging=false;
			GLEventManager::dispatcherList[RIGHT_MOUSE_UP]->dispatch(MouseEvent{x,y});
			break;
		case 30:
			GLEventManager::dispatcherList[MOUSE_WHEEL]->dispatch(MouseEvent{x,y,0,0,1});
			break;
		case 40:
			GLEventManager::dispatcherList[MOUSE_WHEEL]->dispatch(MouseEvent{x,y,0,0,-1});
			break;
	}
}

void MouseDragHandler(int x,int y){

	if(L_draging){
		GLEventManager::dispatcherList[MOUSE_DRAG]->dispatch(MouseEvent{x,y,L_startDragX,L_startDragY});
	}
	if(R_draging){
		GLEventManager::dispatcherList[RIGHT_MOUSE_DRAG]->dispatch(MouseEvent{x,y,R_startDragX,R_startDragY});
	}
	if(M_draging){
		GLEventManager::dispatcherList[MIDDLE_MOUSE_DRAG]->dispatch(MouseEvent{x,y,M_startDragX,M_startDragY});
	}
}

void MouseMoveHandler(int x,int y){

	GLEventManager::dispatcherList[MOUSE_MOVE]->dispatch(MouseEvent{x,y});

}

void KeyDownHandler(unsigned char key, int x, int y){
	GLEventManager::dispatcherList[KEY_DOWN]->dispatch(KeyBoardEvent{key,x,y});
}
void KeyUpHandler(unsigned char key, int x, int y){
	GLEventManager::dispatcherList[KEY_UP]->dispatch(KeyBoardEvent{key,x,y});
}

void SpecialKeyDownHandler(int key, int x, int y){
	unsigned char key_id;
	if(key<100)key_id=key+100;//for F1~F12
	if(key>=100)key_id=key+43;//for control
	GLEventManager::dispatcherList[KEY_DOWN]->dispatch(KeyBoardEvent{key_id,x,y});
}
void SpecialKeyUpHandler(int key, int x, int y){
	unsigned char key_id;
	if(key<100)key_id=key+100;//for F1~F12
	if(key>=100)key_id=key+43;//for control
	GLEventManager::dispatcherList[KEY_UP]->dispatch(KeyBoardEvent{key_id,x,y});
}
