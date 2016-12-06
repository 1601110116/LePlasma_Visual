///*
// * Figure.cpp
// *
// *  Created on: Apr 28, 2016
// *      Author: zlstudio
// */
//
//#include <Figure.h>
//
//Figure::Figure() {
//	// TODO Auto-generated constructor stub
//
//}
//
//Figure::~Figure() {
//	// TODO Auto-generated destructor stub
//}
//
//Figure::Figure(const char* name, Plot& link):content(&link),title(name){
//
//	addChild(link);
//	addChild(box);
//
//	setSize(link.Width+170,link.Height+120,false);
//
//}
//
//void Figure::enterFrame_evt(const Event& evt){
//	content->viewPhi+=AUTO_ROTATE_SPEED;
//}
//
//void Figure::appendDescription(const string& str){
//	if(describe.size()>0){
//		describe+='\n';
//	}
//	describe+=str;
//}
//
//void Figure::clearDescription(){
//	describe.clear();
//}
//void Figure::setSize(int w,int h,bool plotRedraw){
//	Width=w;
//	Height=h;
//
//	if(plotRedraw)content->setSize(Width-170,Height-120);
//
//	graphics.printAt(30+(content->Width)/2,2,title);
//	graphics.printAt(2,Height-120,describe);
//
//	content->x=30;
//	content->y=30;
//
//	graphics.printAt(30+content.Width/2,45+content.Height)
//
//	xlabel.y=;
//	xlabel.x=-xlabel.width/2;
//
//	ylabel.x=10;
//	ylabel.y=30+content.Height/2+ylabel.height/2;
//	ylabel.rotationZ=-90;
//
//}
//
//
