/*
 * Plus.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef PLUS_H_
#define PLUS_H_

#include "Function.h"

#include <vector>

class Plus: public Function {
private:
	Function * f1; Function * f2;
	std::vector<Function*> loose;//pointers to be deleted
public:
	Plus(Function * f1p, Function * f2p);
	~Plus(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	double eval(){return (*f1).eval()+(*f2).eval();}
	Function * diff(Function * var){
		Plus * p = new Plus((*f1).diff(var),(*f2).diff(var));
		loose.push_back(p);
		return p;
	}
};

#endif /* PLUS_H_ */
