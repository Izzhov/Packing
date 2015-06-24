/*
 * Times.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef TIMES_H_
#define TIMES_H_

#include "Function.h"
#include "Plus.h"

#include <vector>

class Times: public Function {
private:
	Function * f1; Function * f2;
	std::vector<Function*> loose;//pointers to be deleted
public:
	Times(Function * f1p, Function * f2p);
	~Times(){for(int i=0; i<loose.size(); i++) delete loose.at(i);}
	double eval(){return (*f1).eval()*(*f2).eval();}
	Function * diff(Function * var){
		Times * t1 = new Times((*f1).diff(var),f2);
		Times * t2 = new Times(f1,(*f2).diff(var));
		Plus * p = new Plus(t1,t2);
		loose.push_back(t1); loose.push_back(t2); loose.push_back(p);
		return p;
	}
};

#endif /* TIMES_H_ */
