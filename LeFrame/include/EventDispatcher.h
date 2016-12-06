/*
 * EventDispatcher.h
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#ifndef EVENTDISPATCHER_H_
#define EVENTDISPATCHER_H_

#include <vector>
#include <Event.h>
#include <EventListener.h>

using namespace std;

class EventDispatcher {
	public:
		EventDispatcher();
		virtual ~EventDispatcher();

		vector<EventListener*> eventQueue;
		vector<void (EventListener::*)(const MouseEvent&)> mouseListenerQueue;
		vector<void (EventListener::*)(const KeyBoardEvent&)> keyBoardListenerQueue;
		vector<void (EventListener::*)(const Event&)> EventListenerQueue;
		vector<void (EventListener::*)(const WindowEvent&)> WindowListenerQueue;

		void dispatch(const MouseEvent&);
		void dispatch(const KeyBoardEvent&);
		void dispatch(const Event&);
		void dispatch(const WindowEvent&);

};

#endif /* EVENTDISPATCHER_H_ */
