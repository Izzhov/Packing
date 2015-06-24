/*
 * Fcon.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef FCON_H_
#define FCON_H_

#include "Variable.h"
#include "Const.h"
#include "Plus.h"
#include "Times.h"
#include "Minus.h"
#include "Div.h"
#include "IntPow.h"
#include "DPow.h"

#include <vector>

class Fcon {
private:
	std::vector<Function*> loose;//pointers to be deleted
public:
	Fcon(){}
	~Fcon(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	Const * c(double con){
		Const * co = new Const(con);
		loose.push_back(co);
		return co;
	}
	Plus * plus(Function * f1, Function * f2){
		Plus * pl = new Plus(f1,f2);
		loose.push_back(pl);
		return pl;
	}
	Times * times(Function *f1, Function * f2){
		Times * ti = new Times(f1,f2);
		loose.push_back(ti);
		return ti;
	}
	Minus * minus(Function * f1, Function * f2){
		Minus * mi = new Minus(f1,f2);
		loose.push_back(mi);
		return mi;
	}
	Div * div(Function * f1, Function * f2){
		Div * di = new Div(f1,f2);
		loose.push_back(di);
		return di;
	}
	IntPow * intpow(Function * f1, int pw){
		IntPow * ip = new IntPow(f1,pw);
		loose.push_back(ip);
		return ip;
	}
	DPow * dpow(Function * f1, double pw){
		DPow * dp = new DPow(f1,pw);
		loose.push_back(dp);
		return dp;
	}
};

#endif /* FCON_H_ */
