/*
 * Minus.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef MINUS_H_
#define MINUS_H_

#include "Function.h"

#include <vector>

class Minus: public Function {
private:
	Function * f1; Function * f2;
	std::vector<Function*> loose;//pointers to be deleted
public:
	Minus(Function * f1p, Function * f2p);
	~Minus(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	double eval(){return (*f1).eval()-(*f2).eval();}
	Function * diff(Function * var){
		Minus * m = new Minus((*f1).diff(var),(*f2).diff(var));
		loose.push_back(m);
		return m;
	}
};

#endif /* MINUS_H_ */
