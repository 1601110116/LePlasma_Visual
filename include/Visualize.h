/*
 * Visualize.h
 *
 *  Created on: Mar 10, 2016
 *      Author: zlstudio
 */

#ifndef VISUALIZE_H_
#define VISUALIZE_H_

#include <LePlasma.h>
#include <fstream>
#include <GridViewer.h>
#include <Plot.h>
#include <Stage.h>
#include <TextField.h>
#include <XYZIndicator.h>
#include <iostream>


using namespace std;

class Plot;
class HostGrid;
class GridViewer;

class Visualize: public Stage {
	public:
		Visualize(Grid*,int ID);
		Visualize(Grid*);
		Visualize(HostGrid*);
		virtual ~Visualize();

		void drawGrid();
		void drawParticles();

		void enterFrame_evt(const Event&);
		void resize_evt(const WindowEvent&);

		void mouseDrag_evt(const MouseEvent& evt);
		void mouseDown_evt(const MouseEvent& evt);
		void mouseUp_evt(const MouseEvent& evt);

		void mouseWheel_evt(const MouseEvent&);
		void keyDown_evt(const KeyBoardEvent&);

		void (*idleFn)();


	private:



		Sprite GridViewPort;

		Plot ExPlot, ExDistribute;

		TextField details;

		GridViewer gridViewer;

		XYZIndicator coordinateIndicator;

		int rx,ry;

		Grid* aimGrid;
		HostGrid* Hostgrid;

		bool isHostView;

		void initAnalysisView();
};

#endif /* VISUALIZE_H_ */
