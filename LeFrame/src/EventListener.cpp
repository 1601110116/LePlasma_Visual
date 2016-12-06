/*
 * EventListener.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#include <EventListener.h>
#include <GLEventManager.h>

EventListener::EventListener() {
	// TODO Auto-generated constructor stub

}

EventListener::~EventListener() {
	// TODO Auto-generated destructor stub
}

void EventListener::addEventListener(eventType mouse_event_type,void (EventListener::*listener)(const MouseEvent&)){
	GLEventManager::registerListener(this,mouse_event_type,listener);
}
void EventListener::addEventListener(eventType key_event_type,void (EventListener::*listener)(const KeyBoardEvent&)){
	GLEventManager::registerListener(this,key_event_type,listener);
}
void EventListener::addEventListener(eventType event_type,void (EventListener::*listener)(const Event&)){
	GLEventManager::registerListener(this,event_type,listener);
}
void EventListener::addEventListener(eventType window_event_type,void (EventListener::*listener)(const WindowEvent&)){
	GLEventManager::registerListener(this,window_event_type,listener);
}

void EventListener::removeEventListener(eventType type){
	GLEventManager::removeListener(this,type);
}

void EventListener::mouseDown_evt(const MouseEvent& evt){}
void EventListener::mouseUp_evt(const MouseEvent& evt){}
void EventListener::mouseDrag_evt(const MouseEvent& evt){}

void EventListener::rightMouseDown_evt(const MouseEvent& evt){}
void EventListener::rightMouseUp_evt(const MouseEvent& evt){}
void EventListener::rightMouseDrag_evt(const MouseEvent& evt){}

void EventListener::middleMouseDown_evt(const MouseEvent& evt){}
void EventListener::middleMouseUp_evt(const MouseEvent& evt){}
void EventListener::middleMouseDrag_evt(const MouseEvent& evt){}

void EventListener::mouseWheel_evt(const MouseEvent& evt){}

void EventListener::mouseMove_evt(const MouseEvent& evt){}

void EventListener::keyDown_evt(const KeyBoardEvent& evt){}

void EventListener::keyUp_evt(const KeyBoardEvent& evt){}

void EventListener::enterFrame_evt(const Event& evt){}

void EventListener::resize_evt(const WindowEvent& evt){}
