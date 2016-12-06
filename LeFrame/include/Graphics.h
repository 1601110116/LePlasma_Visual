/*
 * Graphics.h
 *
 *  Created on: Apr 23, 2016
 *      Author: zlstudio
 */

#ifndef OPENGL_GRAPHICS_H_
#define OPENGL_GRAPHICS_H_

#include <GL/glut.h>

#include <AScompat.h>

#include <Rect.h>

#include <vector>

#include <string>
#include <math.h>

using namespace std;

inline void moveCoordinate(double x,double y,double z){
	glLoadIdentity();
	glTranslated(x,y,z);
}
inline void seeFrom(double r,double theta,double phi){
	gluLookAt(
			r*sin(phi/180*M_PI)*cos(theta/180*M_PI),
			r*sin(theta/180*M_PI),
			r*cos(phi/180*M_PI)*cos(theta/180*M_PI),
			0, 0, 0,
			-sin(phi/180*M_PI)*sin(theta/180*M_PI),
			cos(theta/180*M_PI),
			-cos(phi/180*M_PI)*sin(theta/180*M_PI)
	);
}

inline int getFontHeight(void* font){
	if(font==GLUT_BITMAP_8_BY_13) {
		return 13;
	}else if(font==GLUT_BITMAP_9_BY_15){
		return 15;
	}else if(font==GLUT_BITMAP_HELVETICA_10){
		return 10;
	}else if(font==GLUT_BITMAP_HELVETICA_12){
		return 12;
	}else if(font==GLUT_BITMAP_HELVETICA_18){
		return 18;
	}else if(font==GLUT_BITMAP_TIMES_ROMAN_10){
		return 10;
	}else if(font==GLUT_BITMAP_TIMES_ROMAN_24){
		return 24;
	}
	return 0;
}

struct Stroke{
		double x,y,z;
		bool start;
		bool changeStyle;
};

struct LineStyle{
		int width;
		unsigned char r,g,b,a;
};

struct Text{
		string text;
		void* font;
		int x,y;
		unsigned char r,g,b;
		int h;
};
struct Point{
		double x,y,z;
		bool changeStyle;
};
class Graphics {
	public:
		Graphics();
		virtual ~Graphics();

		///Clear
		void clear();

		///Points
		void pointStyle(int width);
		void pointStyle(int width,unsigned int color,double alpha);
		void drawPoint(double x,double y);
		void drawPoint(double x,double y,double z);


		///2D-Lines
		void lineStyle(int width);
		void lineStyle(int width,unsigned int color,double alpha);
		void moveTo(double x,double y);
		void lineTo(double x,double y);

		///3D-Lines
		void moveTo(double x,double y,double z);
		void lineTo(double x,double y,double z);

		///2D-Drawer
		void drawRect(double x,double y,double w,double h);
		void drawRect(const Rect &rect);

		///3D-Drawer
		void drawBox(double x,double y,double z,double w,double h,double l);
		void drawCube(double x,double y,double z,double a);


		///Text Engine
		/* There are 7 bitmap fonts available in GLUT. They are
		    GLUT_BITMAP_8_BY_13
		    GLUT_BITMAP_9_BY_15
		    GLUT_BITMAP_TIMES_ROMAN_10
		    GLUT_BITMAP_TIMES_ROMAN_24
		    GLUT_BITMAP_HELVETICA_10
		    GLUT_BITMAP_HELVETICA_12
		    GLUT_BITMAP_HELVETICA_18 */
		void textStyle(Color color);
		void textStyle(void* font,Color color);
		void printAt(int x,int y,const char* text);
		void printAt(int x,int y,const string &s);

	private:
		unsigned char TextColor[3];
		int Texth;
		void* TextFont;

		vector<Point> points;
		vector<Stroke> strokes;
		vector<LineStyle> lineStyles;
		vector<LineStyle> pointStyles;

		vector<Text> texts;



	protected:
		void draw();

		friend class DisplayObject;
		friend class TextField;


};

#endif /* OPENGL_GRAPHICS_H_ */
