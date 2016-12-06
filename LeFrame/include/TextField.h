/*
 * TextField.h
 *
 *  Created on: Apr 28, 2016
 *      Author: zlstudio
 */

#ifndef TEXTFIELD_H_
#define TEXTFIELD_H_

#include <DisplayObject.h>
#include <string>

using namespace std;

class TextField :public DisplayObject{
	public:
		TextField();
		virtual ~TextField();

		void setTextStyle(Color c);
		void setTextStyle(void* font, Color c);

		void setText(const string &str);
		string& getText();


		int textWidth();
		int textHeight();

	private:
		int Width,Height,rows;
		string text;

};

#endif /* TEXTFIELD_H_ */
