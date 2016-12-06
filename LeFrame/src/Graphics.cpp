/*
 * Graphics.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: zlstudio
 */

#include <Graphics.h>

inline void setColor(Color color){
	unsigned char r,g,b;
	phraseColor(color,r,g,b);
	glColor3ub(r,g,b);
}

Graphics::Graphics():TextFont(GLUT_BITMAP_HELVETICA_12),Texth(12){
	// TODO Auto-generated constructor stub

}

Graphics::~Graphics() {
	// TODO Auto-generated destructor stub
}

///Drawers
bool changeNext=false;
bool P_changeNext=false;

void Graphics::lineStyle(int width,unsigned int color,double alpha){


	unsigned char r,g,b;
	unsigned char a=alpha*255;
	phraseColor(color,r,g,b);

	if(changeNext){
		//means last setting is not used
		int last=lineStyles.size()-1;
		lineStyles[last].width=width;
		lineStyles[last].r=r;
		lineStyles[last].g=g;
		lineStyles[last].b=b;
		lineStyles[last].a=a;

	}else{
		lineStyles.push_back(LineStyle{width,r,g,b,a});
	}

	changeNext=true;
}

void Graphics::pointStyle(int width,unsigned int color,double alpha){


	unsigned char r,g,b;
	unsigned char a=alpha*255;
	phraseColor(color,r,g,b);

	if(P_changeNext){
		//means last setting is not used
		int last=pointStyles.size()-1;
		pointStyles[last].width=width;
		pointStyles[last].r=r;
		pointStyles[last].g=g;
		pointStyles[last].b=b;
		pointStyles[last].a=a;

	}else{
		pointStyles.push_back(LineStyle{width,r,g,b,a});
	}

	P_changeNext=true;
}

void Graphics::clear(){
	//	stroke includes no pointer. So no delete is needed
	strokes.clear();
	lineStyles.clear();
	pointStyles.clear();
	points.clear();
	texts.clear();
	P_changeNext=false;
	changeNext=false;
}

void Graphics::drawRect(double x,double y,double width,double height){

	strokes.push_back(Stroke{x,-y,0,true,changeNext});
	strokes.push_back(Stroke{x+width,-y,0});
	strokes.push_back(Stroke{x+width,-y-height,0});
	strokes.push_back(Stroke{x,-y-height,0});
	strokes.push_back(Stroke{x,-y,0});

	changeNext=false;

}

void Graphics::drawRect(const Rect &rect){
	drawRect(rect.x,rect.y,rect.width,rect.width);
}

void Graphics::drawPoint(double x,double y){

	points.push_back(Point{x,-y,0,P_changeNext});
	P_changeNext=false;

}

void Graphics::drawPoint(double x,double y,double z){

	points.push_back(Point{x,-y,-z,P_changeNext});
	P_changeNext=false;

}

void Graphics::moveTo(double x,double y){

	strokes.push_back(Stroke{x,-y,0,true,changeNext});
	changeNext=false;

}

void Graphics::lineTo(double x,double y){

	strokes.push_back(Stroke{x,-y,0,false,changeNext});
	changeNext=false;

}

void Graphics::moveTo(double x,double y,double z){

	strokes.push_back(Stroke{x,-y,-z,true,changeNext});
	changeNext=false;

}

void Graphics::lineTo(double x,double y,double z){

	strokes.push_back(Stroke{x,-y,-z,false,changeNext});
	changeNext=false;

}

///3D-Drawers
void Graphics::drawBox(double x,double y,double z,double width,double height,double length){
	strokes.push_back(Stroke{x,-y,-z,true,changeNext});
	strokes.push_back(Stroke{x+width,-y,-z});
	strokes.push_back(Stroke{x+width,-y-height,-z});
	strokes.push_back(Stroke{x,-y-height,-z});
	strokes.push_back(Stroke{x,-y,-z});

	strokes.push_back(Stroke{x,-y,-z-length});
	strokes.push_back(Stroke{x+width,-y,-z-length});
	strokes.push_back(Stroke{x+width,-y-height,-z-length});
	strokes.push_back(Stroke{x,-y-height,-z-length});
	strokes.push_back(Stroke{x,-y,-z-length});


	strokes.push_back(Stroke{x+width,-y,-z,true});
	strokes.push_back(Stroke{x+width,-y,-z-length});

	strokes.push_back(Stroke{x+width,-y-height,-z,true});
	strokes.push_back(Stroke{x+width,-y-height,-z-length});

	strokes.push_back(Stroke{x,-y-height,-z,true});
	strokes.push_back(Stroke{x,-y-height,-z-length});

	changeNext=false;
}

void Graphics::drawCube(double x,double y,double z,double a){
	drawBox(x,y,z,a,a,a);
}


////Temporary Text Agent

void Graphics::textStyle(Color color){
	phraseColor(color,TextColor[0],TextColor[1],TextColor[2]);
}
void Graphics::textStyle(void* font,Color color){
	TextFont=font;
	phraseColor(color,TextColor[0],TextColor[1],TextColor[2]);

	Texth=getFontHeight(font);

}
void Graphics::printAt(int x,int y,const char* text){

	string s(text);

	texts.push_back(Text{s,TextFont,x,-y,TextColor[0],TextColor[1],TextColor[2],Texth});

}

void Graphics::printAt(int x,int y,const string &s){
	texts.push_back(Text{s,TextFont,x,-y,TextColor[0],TextColor[1],TextColor[2],Texth});
}


///Central Drawer
void Graphics::draw(){

	//draw lines
	int styleIndex=0;
	glColor4ub(0,0,0,0);
	if(strokes.size()>0){

		glBegin(GL_LINE_STRIP);
		for (unsigned int i = 0; i < strokes.size(); ++i) {
			if(strokes[i].start&&i!=0){
				glEnd();
				glBegin(GL_LINE_STRIP);
			}
			if(strokes[i].changeStyle){
				glColor4ub(
						lineStyles[styleIndex].r,
						lineStyles[styleIndex].g,
						lineStyles[styleIndex].b,
						lineStyles[styleIndex].a
				);
				styleIndex++;
			}
			glVertex3d(strokes[i].x,strokes[i].y,strokes[i].z);
		}
		glEnd();
	}

	if(points.size()>0){
		styleIndex=0;
		glBegin(GL_POINTS);
		for (unsigned int i = 0; i < points.size(); ++i) {
			if(points[i].changeStyle){
				glColor4ub(
						pointStyles[styleIndex].r,
						pointStyles[styleIndex].g,
						pointStyles[styleIndex].b,
						pointStyles[styleIndex].a
				);
				styleIndex++;
			}
			glVertex3d(points[i].x,points[i].y,points[i].z);
		}
		glEnd();
	}

	//draw text
	for (unsigned int i = 0; i < texts.size(); ++i) {
		glColor3ub(
				texts[i].r,
				texts[i].g,
				texts[i].b
		);
		int startH=texts[i].y-texts[i].h;
		glRasterPos2f(texts[i].x,startH);
		for (unsigned int j = 0; j < texts[i].text.size(); ++j) {
			if(texts[i].text[j]=='\n'){
				startH-=texts[i].h;
				glRasterPos2f(texts[i].x,startH);
			}else{
				glutBitmapCharacter(texts[i].font, texts[i].text[j]);
			}
		}
	}

}
