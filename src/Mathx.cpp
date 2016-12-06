/*
 * Mathx.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: zlstudio
 */

#include <Mathx.h>


Mathx::Mathx() {
	// TODO Auto-generated constructor stub

}

Mathx::~Mathx() {
	// TODO Auto-generated destructor stub
}

Arrange Mathx::Tick125(double s,double e){
	double base10=pow(10,floor(log10(e-s)));
	double tick1=base10/10.0;
	double tick2=base10/5.0;
	double tick5=base10/2.0;

	double baseStart=floorFloat(s,-floor(log10(e-s)));
	double baseEND=ceilFloat(e,-floor(log10(e-s)));


	double realStart1=baseStart+floor((s-baseStart)/tick1)*tick1;
	double realEnd1=baseEND-floor((baseEND-e)/tick1)*tick1;
	int ticknumbers1=(realEnd1-realStart1)/tick1;

	double realStart2=baseStart+floor((s-baseStart)/tick2)*tick2;
	double realEnd2=baseEND-floor((baseEND-e)/tick2)*tick2;
	int ticknumbers2=(realEnd2-realStart2)/tick2;

	double realStart5=baseStart+floor((s-baseStart)/tick5)*tick5;
	double realEnd5=baseEND-floor((baseEND-e)/tick5)*tick5;
	int ticknumbers5=(realEnd5-realStart5)/tick5;

	int ticknumbers10=(baseEND-baseStart)/base10;


	int dev1=abs(10-ticknumbers1);
	int dev2=abs(10-ticknumbers2);
	int dev5=abs(10-ticknumbers5);
	int dev10=abs(10-ticknumbers10);

	if(dev1<dev2&&dev1<dev10&&dev1<dev5&&ticknumbers1>4)return Arrange {realStart1,realEnd1,tick1,base10};
	if(dev2<dev5&&dev2<dev10&&ticknumbers2>4)return Arrange {realStart2,realEnd2,tick2,base10};
	if(dev5<dev10&&ticknumbers5>4)return Arrange {realStart5,realEnd5,tick5,base10};

	return Arrange {baseStart,baseEND,base10,base10};
}

double Mathx::ceilFloat(double n,int d)
{
	double e=pow(10,d);
	return ceil(n*e)/e;
}
double Mathx::floorFloat(double n,int d)
{
	double e=pow(10,d);
	return floor(n*e)/e;
}
string Mathx::showFloat(double n,int d){

	ostringstream tmpchar;

	tmpchar<<n;

	return tmpchar.str();
}

string Mathx::DoubleToString(double d)
{
	string str;
	stringstream ss;
	ss<<d;
	ss>>str;
	return str;
}
