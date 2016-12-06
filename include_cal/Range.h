/*
 * Range.h
 *
 *  Created on: Mar 13, 2016
 *      Author: zlstudio
 */

#ifndef RANGE_H_
#define RANGE_H_

#include<string>
#include<sstream>

using namespace std;

class Range {
	public:
		Range();
		Range(int x1,int x2,int y1,int y2,int z1,int z2);
		virtual ~Range();

		int sx,sy,sz;
		int ex,ey,ez;

		int rangeX() const;
		int rangeY() const;
		int rangeZ() const;

		int Volumn() const;

		Range &operator =(const Range &a);

		string toString();

	private:
		int lx,ly,lz;
		int v;
};

#endif /* RANGE_H_ */
