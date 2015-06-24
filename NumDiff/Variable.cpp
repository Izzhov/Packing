/*
 * Variable.cpp
 *
 *  Created on: Jun 17, 2015
 *      Author: izzhov
 */

#include "Variable.h"

Variable::Variable(double var) {
	set_var(var);
}

void Variable::set_var(double var){
	this->var = var;
}

