/*
 * UIContainer.h
 *
 *  Created on: Apr 24, 2016
 *      Author: zlstudio
 */

#ifndef LEFRAME_UICONTAINER_H_
#define LEFRAME_UICONTAINER_H_


#include <Graphics.h>
#include <Rect.h>
#include <vector>

class Stage;

using namespace std;

class DisplayObject :public Rect{
	public:
		DisplayObject();
		virtual ~DisplayObject();

		void addChild(DisplayObject*);
		void addChild(DisplayObject&);
		void removeChild(DisplayObject*);
		void removeChildAt(int);

		double viewDistance;
		double viewTheta;
		double viewPhi;

		Stage *stage;

		Graphics graphics;

		bool hitTest(double x,double y);

	protected:
		vector<DisplayObject*> childList;
		void redraw(double,double);

		double globalX,globalY;


};



#endif /* LEFRAME_UICONTAINER_H_ */
