/*
 * SimplexFltr.h
 *
 *  Created on: Jul 16, 2015
 *      Author: izzhov
 */

#ifndef SIMPLEXFLTR_H_
#define SIMPLEXFLTR_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Torus.h"
#include "Sphere.h"
#include "Spherocyl.h"
#include "mm.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <cstddef>
#include <vector>

//this class is meant to work in conjunction with ConNetFltr
//assumes number of contacts is > number of d.o.f.'s
//see 7-16-15 notes for detailed explanation of this algorithm
class SimplexFltr {
private:
	std::vector<gsl_vector*> ftis;
	std::vector<gsl_vector*> gjs;
	std::vector<int> bi;
	int c; //number of contacts
	int d;//number of d.o.f.'s

	double signtol;//if abs(num) < signtol, then num is considered equal to zero

	int maxsteps;//maximum number of simplex algo steps you're allowed to do
public:
	//stores fis in class (as ftis) and initializes gjs and bi entries to zeroes
	//also stores c and d
	SimplexFltr(std::vector<gsl_vector*> fis, std::vector<gsl_vector*> flocs,
			int d, int sym);

	//only call this once per class instance!
	//after this finishes, it'll free every gsl_vector in the class
	bool is_floater();
	//checks if applying force sign*e_(unitj) is ok
	bool is_floater(int unitj, int sign);

	//the first column with a negative value in bottom row (for entr-ing var)
	int entering_col();

	//pivots about j-th row and i-th column (makes column a basis vector)
	//precondition: gjs_j^i>signtol
	void pivot(int jj, int ii);

	//prints the tableau in its current state
	void print_tab();
};

#endif /* SIMPLEXFLTR_H_ */
