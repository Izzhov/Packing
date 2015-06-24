/*
 * DPow.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef DPOW_H_
#define DPOW_H_

#include "Function.h"
#include "Const.h"
#include "Times.h"

#include <vector>
#include <cmath>

class DPow: public Function {
private:
	Function * f1; double pw;
	std::vector<Function*> loose;//pointers to be deleted
public:
	DPow(Function * f1p, double power);
	~DPow(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	double eval(){return std::pow((*f1).eval(),pw);}
	Function * diff(Function * var){
		Const * cpw = new Const(pw);
		loose.push_back(cpw);
		DPow * newpow = new DPow(f1,pw-1);
		loose.push_back(newpow);
		Times * t = new Times(cpw,newpow);
		loose.push_back(t);
		Times * t2 = new Times(t,(*f1).diff(var));
		loose.push_back(t2);
		return t2;
	}
};

#endif /* DPOW_H_ */
