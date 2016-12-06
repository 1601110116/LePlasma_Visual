/*
 * AScompat.h
 *
 *  Created on: Apr 23, 2016
 *      Author: zlstudio
 */

#ifndef INCLUDE_ASCOMPAT_H_
#define INCLUDE_ASCOMPAT_H_


#define Max(a,b) a>b?a:b
#define Min(a,b) a<b?a:b

#define DEFAULT_STAGE_WIDTH 500
#define DEFAULT_STAGE_HEIGHT 400

typedef unsigned int Color;

inline void phraseColor(Color color, unsigned char &r,unsigned char &g,unsigned char &b){
	r=(color&0xff0000)>>16;
	g=(color&0x00ff00)>>8;
	b=(color&0x0000ff);
}

#endif /* INCLUDE_ASCOMPAT_H_ */
