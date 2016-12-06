/*
 * Key.h
 *
 *  Created on: Apr 25, 2016
 *      Author: zlstudio
 */

#ifndef LEFRAME_INCLUDE_KEYBOARD_H_
#define LEFRAME_INCLUDE_KEYBOARD_H_

///all inputable charactors are within 32~127.
///OpenGL do not tell the different between number pad and number key, or the function of number pad

//some useful key within 0~31:
#define BACKSPACE 8
#define TAB 9
#define ENTER 13
#define ESCAPE 27
#define SPACE 32
#define DELETE 46 //127

//Then we should map function keys greater than 127.
//here we use 128~256 to keep using unsigned char

///F1~F12 in openGL is S1~S12, now +127
#define F1 128
#define F2 129
#define F3 130
#define F4 131
#define F5 132
#define F6 133
#define F7 134
#define F8 135
#define F9 136
#define F10 137
#define F11 138
#define F12 139
#define F13 140
#define F14 141
#define F15 142

//S100~S109: +43
#define LEFT 143
#define UP 144
#define RIGHT 145
#define DOWN 146
#define PAGE_UP 147
#define PAGE_DOWN 148
#define HOME 149
#define END 150
#define INSERT 151
#define NUMBER_LOCK 152

//S112~S117: +43
#define L_SHIFT 155
#define R_SHIFT 156
#define L_CTRL 157
#define R_CTRL 158
#define L_ALT 159
#define R_ALT 160




#endif /* LEFRAME_INCLUDE_KEYBOARD_H_ */
