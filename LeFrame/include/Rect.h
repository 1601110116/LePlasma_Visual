/*
 * Rect.h
 *
 *  Created on: Apr 21, 2016
 *      Author: zlstudio
 */

#ifndef OPENGL_RECT_H_
#define OPENGL_RECT_H_


class Rect {
	public:
		Rect(double x,double y,double w,double h);
		virtual ~Rect();

		double x,y,width,height;
};

#endif /* OPENGL_RECT_H_ */
