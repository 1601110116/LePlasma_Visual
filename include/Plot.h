/*
 * Plot.h
 *
 *  Created on: Apr 23, 2016
 *      Author: zlstudio
 */

#ifndef SRC_PLOT_H_
#define SRC_PLOT_H_

#define ITEM_INDICATOR 1
#define ITEM_INDICATORX 2
#define ITEM_INDICATORY 3
#define ITEM_POINTS 4
#define ITEM_FUNCTION 5
#define ITEM_TEXT 6

#define FONT_COLOR 0xcccccc
#define AXIS_COLOR 0x0000ff
#define FRAME_COLOR 0x9999ff
#define GRAPH_COLOR 0xffff00

#include <vector>
#include <Mathx.h>
#include <Vector3D.h>
#include <Sprite.h>
#include <AScompat.h>

struct PlotRange{
		double sx,ex,sy,ey;
};

using namespace std;

class Plot :public Sprite{
	public:
		Plot();
		Plot(int w,int h,const char* xlab ,const char* ylab);
		Plot(int w,int h,double sx,double ex,double sy,double ey,const char* xlab ,const char* ylab);
		Plot(int w,int h,double sx,double ex,double sy,double ey,int xlog,int ylog,const char* xlab ,const char* ylab,double tickx,double ticky);
		virtual ~Plot();

		int Width,Height;
		double SX,SY,EX,EY;

		string xlabel;
		string ylabel;

		string xname;
		string yname;

		bool autoRange;

		void setCoordinate(const PlotRange&);

		void setSize(int w,int h);

		void trackData1D(const vector<double>&);

		void drawFunction(double (*fn)(double x), Color c);


	private:
		DisplayObject graph;

		double xscale;
		double yscale;

		int baseX;
		int baseY;

		double basescaleX;
		double basescaleY;

		double xtick,ytick;

		bool noRange;

		//Range: L,R,T,B
		PlotRange range;

};



#endif /* SRC_PLOT_H_ */
