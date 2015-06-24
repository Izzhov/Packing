/*
 * Variable.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef VARIABLE_H_
#define VARIABLE_H_

#include "Function.h"
#include "Const.h"

#include <boost/shared_ptr.hpp>

#include <vector>

class Variable: public Function {
private:
	double var;
	std::vector<Function*> loose;//pointers to be deleted
public:
	Variable(){}
	Variable(double var);
	~Variable(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	void set_var(double var);
	double eval(){return var;}
	Function * diff(Function * var){
		if(var==this){
			Const * v = new Const(1);
			loose.push_back(v);
			return v;
		}
		else{
			Const * v = new Const(0);
			loose.push_back(v);
			return v;
		}
	}
};

#endif /* VARIABLE_H_ */
