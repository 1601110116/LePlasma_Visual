/*
 * TextField.cpp
 *
 *  Created on: Apr 28, 2016
 *      Author: zlstudio
 */

#include <GL/freeglut_std.h>
#include <TextField.h>
#include <string>

TextField::TextField() {
	// TODO Auto-generated constructor stub
	graphics.textStyle(GLUT_BITMAP_HELVETICA_12,0x000000);
	graphics.printAt(0,0,string(""));
}

TextField::~TextField() {
	// TODO Auto-generated destructor stub
}

int TextField::textHeight(){
	return Height;
}
int TextField::textWidth(){

}
void TextField::setTextStyle(void* font,Color c){
	unsigned char r,g,b;
	phraseColor(c,r,g,b);
	graphics.texts[0].font=font;
	graphics.texts[0].r=r;
	graphics.texts[0].g=g;
	graphics.texts[0].b=b;
	graphics.texts[0].h=getFontHeight(font);

	graphics.textStyle(font,c);

	Height=graphics.Texth*rows;
}

void TextField::setTextStyle(Color c){
	unsigned char r,g,b;
	phraseColor(c,r,g,b);
	graphics.texts[0].r=r;
	graphics.texts[0].g=g;
	graphics.texts[0].b=b;

	graphics.textStyle(c);
}

void TextField::setText(const string &str){
	text=str;
	graphics.texts[0].text=str;
	rows=1;

	for (int i = 0; i < str.size(); ++i) {
		if(str[i]=='\n')rows++;
	}

	Height=graphics.Texth*rows;
}

string& TextField::getText(){
	return text;
}
