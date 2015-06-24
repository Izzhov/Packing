/*
 * IntPow.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef INTPOW_H_
#define INTPOW_H_

#include "Function.h"
#include "Const.h"
#include "Times.h"

#include <vector>
#include <cmath>

class IntPow: public Function {
private:
	Function * f1; int pw;
	std::vector<Function*> loose;//pointers to be deleted
public:
	IntPow(Function * f1p, int power);
	~IntPow(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	double eval(){return std::pow((*f1).eval(),pw);}
	Function * diff(Function * var){
		if(pw==2){
			Const * cpw = new Const(pw);
			loose.push_back(cpw);
			Times * t = new Times(cpw,f1);
			loose.push_back(t);
			Times * t2 = new Times(t,(*f1).diff(var));
			loose.push_back(t2);
			return t2;
		}
		else{
			Const * cpw = new Const(pw);
			loose.push_back(cpw);
			IntPow * newpow = new IntPow(f1,pw-1);
			loose.push_back(newpow);
			Times * t = new Times(cpw,newpow);
			loose.push_back(t);
			Times * t2 = new Times(t,(*f1).diff(var));
			loose.push_back(t2);
			return t2;
		}
	}
};

#endif /* INTPOW_H_ */
