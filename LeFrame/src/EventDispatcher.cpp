/*
 * EventDispatcher.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#include <EventDispatcher.h>

EventDispatcher::EventDispatcher() {
	// TODO Auto-generated constructor stub

}

EventDispatcher::~EventDispatcher() {
	// TODO Auto-generated destructor stub
}

void EventDispatcher::dispatch(const MouseEvent& event){
	for (unsigned int i = 0; i < mouseListenerQueue.size(); ++i) {
		(eventQueue[i]->*mouseListenerQueue[i])(event);
	}
}
void EventDispatcher::dispatch(const KeyBoardEvent& event){
	for (unsigned int i = 0; i < keyBoardListenerQueue.size(); ++i) {
		(eventQueue[i]->*keyBoardListenerQueue[i])(event);
	}
}

void EventDispatcher::dispatch(const Event& event){
	for (unsigned int i = 0; i < EventListenerQueue.size(); ++i) {
		(eventQueue[i]->*EventListenerQueue[i])(event);
	}
}
void EventDispatcher::dispatch(const WindowEvent& event){
	for (unsigned int i = 0; i < WindowListenerQueue.size(); ++i) {
		(eventQueue[i]->*WindowListenerQueue[i])(event);
	}
}
