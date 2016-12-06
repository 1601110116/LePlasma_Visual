///*
// * Figure.h
// *
// *  Created on: Apr 28, 2016
// *      Author: zlstudio
// */
//
//#ifndef SRC_FIGURE_H_
//#define SRC_FIGURE_H_
//
//#include <Plot.h>
//#include <Sprite.h>
//#include <string>
//#include <Event.h>
//
//#define AUTO_ROTATE_SPEED 0.0175
//
//using namespace std;
//
//class Figure :public Sprite{
//	public:
//		Figure();
//		Figure(const char* name, Plot&);
//		virtual ~Figure();
//
//		static bool AutoRotate;
//
//		int Height,Width;
//
//		Plot* content;
//
//		void appendDescription(const string& str);
//		void clearDescription();
//		void setSize(int w,int h,bool plotRedraw);
//
//	private:
//		string title;
//		DisplayObject box;
//		string describe;
//
//		string* xlabel;
//		string* ylabel;
//
//		int at;
//
//		void enterFrame_evt(const Event&);
//};
//
//#endif /* SRC_FIGURE_H_ */
