/*
 * EventListener.h
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#ifndef EVENTLISTENER_H_
#define EVENTLISTENER_H_

#include <Event.h>

class EventListener {
	public:
		EventListener();
		virtual ~EventListener();

		void addEventListener(eventType,void (EventListener::*listener)(const MouseEvent&));
		void addEventListener(eventType,void (EventListener::*listener)(const KeyBoardEvent&));
		void addEventListener(eventType,void (EventListener::*listener)(const Event&));
		void addEventListener(eventType,void (EventListener::*listener)(const WindowEvent&));



		void removeEventListener(eventType);


		virtual void mouseDown_evt(const MouseEvent& evt);
		virtual void mouseUp_evt(const MouseEvent& evt);
		virtual void mouseDrag_evt(const MouseEvent& evt);

		virtual void rightMouseDown_evt(const MouseEvent& evt);
		virtual void rightMouseUp_evt(const MouseEvent& evt);
		virtual void rightMouseDrag_evt(const MouseEvent& evt);

		virtual void middleMouseDown_evt(const MouseEvent& evt);
		virtual void middleMouseUp_evt(const MouseEvent& evt);
		virtual void middleMouseDrag_evt(const MouseEvent& evt);

		virtual void mouseWheel_evt(const MouseEvent& evt);

		virtual void mouseMove_evt(const MouseEvent& evt);

		virtual void keyDown_evt(const KeyBoardEvent& evt);

		virtual void keyUp_evt(const KeyBoardEvent& evt);

		virtual void enterFrame_evt(const Event& evt);

		virtual void resize_evt(const WindowEvent& evt);
};

#endif /* EVENTLISTENER_H_ */
