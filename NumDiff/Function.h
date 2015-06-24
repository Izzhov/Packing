/*
 * Function.h
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#ifndef FUNCTION_H_
#define FUNCTION_H_

#include <boost/shared_ptr.hpp>

class Function {
public:
	virtual ~Function(){}
	virtual double eval() = 0;
	virtual Function * diff(Function * var) = 0;
};

#endif /* FUNCTION_H_ */
