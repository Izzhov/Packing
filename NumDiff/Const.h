/*
 * Const.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef CONST_H_
#define CONST_H_

#include "Function.h"

#include <boost/shared_ptr.hpp>
#include <vector>

class Const: public Function {
private:
	double val;
	std::vector<Function*> loose;//pointers to be deleted
public:
	Const(double val);
	~Const(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	double eval(){return val;}
	Function * diff(Function * var){
		Const * v = new Const(0);
		loose.push_back(v);
		return v;
	}
};

#endif /* CONST_H_ */
