/*
 * CompareHPs.h
 *
 *  Created on: Jul 10, 2015
 *      Author: izzhov
 */

#ifndef COMPAREHPS_H_
#define COMPAREHPS_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <sstream>

#include <mm.h>
#include <Torus.h>
#include <SpherocylDOF.h>
#include <HarmPot.h>
#include <HarmPotNbrList.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <vector>
#include <cmath>

class CompareHPs {
private:
	double alph; double alph2;
	double eps; double eps2;
	double sig;
	int imax; //see conjgrad paper for meaning of these
	HarmPot<Torus<SpherocylDOF> >& pot; Torus<SpherocylDOF>& box;
	HarmPotNbrList<Torus<SpherocylDOF>,SpherocylDOF>& pot2;
	Torus<SpherocylDOF>& box2;
public:
	CompareHPs(HarmPot<Torus<SpherocylDOF> >& ppot, Torus<SpherocylDOF>& bbox,
			HarmPotNbrList<Torus<SpherocylDOF>,SpherocylDOF>& ppot2, Torus<SpherocylDOF>& bbox2);
	bool minimize();

	void move(double alphsig, int num);
	void translate(double alphsig, int num);
	void rotate(double alphsig, int num);

	void store_deriv(gsl_vector * r, int scale, int num);//scale to store pos or neg
	void load_deriv(gsl_vector * r, int scale, int num);
	//returns the vector^2 but with every sym-th element mult'd by L
	//(to fix the units on the forces - see7-9-15 whiteback)
	double unit_fix(gsl_vector * r, int num);
};

#endif /* COMPAREHPS_H_ */
