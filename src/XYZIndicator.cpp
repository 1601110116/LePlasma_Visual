/*
 * XYZIndicator.cpp
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#include <XYZIndicator.h>

#define scale 20

XYZIndicator::XYZIndicator() {

	graphics.lineStyle(0,0xff0000,1);
	graphics.moveTo(0,0,0);
	graphics.lineTo(scale,0,0);

	graphics.lineStyle(0,0x00ff00,1);

	graphics.moveTo(0,0,0);
	graphics.lineTo(0,scale,0);

	graphics.lineStyle(0,0x0000ff,1);
	graphics.moveTo(0,0,0);
	graphics.lineTo(0,0,scale);

}

XYZIndicator::~XYZIndicator() {
	// TODO Auto-generated destructor stub
}

