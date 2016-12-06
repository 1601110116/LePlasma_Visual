/*
 * GLEventManager.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#include <GLEventManager.h>

vector<EventDispatcher*> GLEventManager::dispatcherList;

GLEventManager::GLEventManager() {
	for (int i = 0; i < EVENT_TYPES; ++i) {
		dispatcherList.push_back(new EventDispatcher());
	}
}

GLEventManager::~GLEventManager() {
	// TODO Auto-generated destructor stub
}

void GLEventManager::registerListener(EventListener* target,eventType id,void (EventListener::*listener)(const MouseEvent&)){
	dispatcherList[id]->eventQueue.push_back(target);
	dispatcherList[id]->mouseListenerQueue.push_back(listener);
}

void GLEventManager::registerListener(EventListener* target,eventType id,void (EventListener::*listener)(const KeyBoardEvent&)){
	dispatcherList[id]->eventQueue.push_back(target);
	dispatcherList[id]->keyBoardListenerQueue.push_back(listener);
}
void GLEventManager::registerListener(EventListener* target,eventType id,void (EventListener::*listener)(const Event&)){
	dispatcherList[id]->eventQueue.push_back(target);
	dispatcherList[id]->EventListenerQueue.push_back(listener);
}
void GLEventManager::registerListener(EventListener* target,eventType id,void (EventListener::*listener)(const WindowEvent&)){
	dispatcherList[id]->eventQueue.push_back(target);
	dispatcherList[id]->WindowListenerQueue.push_back(listener);
}

void GLEventManager::removeListener(EventListener* target,eventType id){
	for (unsigned int i = 0; i < dispatcherList[id]->eventQueue.size(); ++i) {
		if(dispatcherList[id]->eventQueue[i]==target){
			dispatcherList[id]->eventQueue.erase(dispatcherList[id]->eventQueue.begin()+i);
			if(dispatcherList[id]->mouseListenerQueue.size()>0){
				dispatcherList[id]->mouseListenerQueue.erase(dispatcherList[id]->mouseListenerQueue.begin()+i);
			}else if(dispatcherList[id]->keyBoardListenerQueue.size()>0){
				dispatcherList[id]->keyBoardListenerQueue.erase(dispatcherList[id]->keyBoardListenerQueue.begin()+i);
			}else if(dispatcherList[id]->EventListenerQueue.size()>0){
				dispatcherList[id]->EventListenerQueue.erase(dispatcherList[id]->EventListenerQueue.begin()+i);
			}else if(dispatcherList[id]->WindowListenerQueue.size()>0){
				dispatcherList[id]->WindowListenerQueue.erase(dispatcherList[id]->WindowListenerQueue.begin()+i);
			}
			return;
		}
	}
}
