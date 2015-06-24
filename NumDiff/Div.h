/*
 * Div.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef DIV_H_
#define DIV_H_

#include "Function.h"
#include "Times.h"
#include "Minus.h"
#include "IntPow.h"

class Div: public Function {
private:
	Function * f1; Function * f2;
	std::vector<Function*> loose;//pointers to be deleted
public:
	Div(Function * f1p, Function * f2p);
	~Div(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	double eval(){return (*f1).eval()/(*f2).eval();}
	Function * diff(Function * var){
		Times * t1 = new Times((*f1).diff(var),f2);
		loose.push_back(t1);
		Times * t2 = new Times(f1,(*f2).diff(var));
		loose.push_back(t2);
		Minus * m = new Minus(t1,t2);
		loose.push_back(m);
		IntPow * ip = new IntPow(f2,2);
		loose.push_back(ip);
		Div * di = new Div(m,ip);
		loose.push_back(di);
		return di;
	}
};

#endif /* DIV_H_ */
