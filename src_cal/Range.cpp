/*
 * Range.cpp
 *
 *  Created on: Mar 13, 2016
 *      Author: zlstudio
 */

#include "Range.h"

Range::Range():
sx(0),ex(0),
sy(0),ey(0),
sz(0),ez(0),
lx(0),ly(0),lz(0),
v(0){}

Range::~Range() {
	// TODO Auto-generated destructor stub
}

Range::Range(int x1,int x2,int y1,int y2,int z1,int z2):
sx(x1),ex(x2),
sy(y1),ey(y2),
sz(z1),ez(z2),
lx(ex-sx),ly(ey-sy),lz(ez-sz),
v(lx*ly*lz)
{}

int Range::rangeX()const{return lx;}
int Range::rangeY()const{return ly;}
int Range::rangeZ()const{return lz;}

int Range::Volumn()const{return v;}

Range &Range::operator =(const Range &a){
	sx=a.sx;ex=a.ex;
	sy=a.sy;ey=a.ey;
	sz=a.sz;ez=a.ez;

	lx=ex-sx;
	ly=ey-sy;
	lz=ez-sz;

	v=lx*ly*lz;
	return *this;
}

string Range::toString(){
	ostringstream tmpchar;

	tmpchar<<"("<<sx<<" -> "<<ex<<", "<<sy<<" -> "<<ey<<", "<<sz<<" -> "<<ez<<")";

	return tmpchar.str();
}

