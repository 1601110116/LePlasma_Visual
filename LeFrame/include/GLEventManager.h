/*
 * GLEventManager.h
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#ifndef GLEVENTMANAGER_H_
#define GLEVENTMANAGER_H_


#include <Event.h>
#include <EventDispatcher.h>
#include <EventListener.h>
//static

using namespace std;

class GLEventManager {
	public:
		GLEventManager();
		virtual ~GLEventManager();

		static vector<EventDispatcher*> dispatcherList;
		static void registerListener(EventListener*,eventType,void (EventListener::*listener)(const MouseEvent&));
		static void registerListener(EventListener*,eventType,void (EventListener::*listener)(const KeyBoardEvent&));
		static void registerListener(EventListener*,eventType,void (EventListener::*listener)(const Event&));
		static void registerListener(EventListener*,eventType,void (EventListener::*listener)(const WindowEvent&));


		static void removeListener(EventListener*,eventType);

};

#endif /* GLEVENTMANAGER_H_ */
