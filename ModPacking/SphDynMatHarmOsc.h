/*
 * SphDynMatHarmOsc.h
 *
 *  Created on: Jun 18, 2015
 *      Author: izzhov
 */

#ifndef SPHDYNMATHARMOSC_H_
#define SPHDYNMATHARMOSC_H_

#include "Torus.h"
#include "Sphere.h"
#include "HarmPot.h"
#include "SphereConNet.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <vector>
#include <cmath>

#include <Fcon.h>
#include <Variable.h>
#include <Function.h>

class SphDynMatHarmOsc {
private:
	Torus<Sphere>& box;
	HarmPot<Torus<Sphere> >& pot;
	SphereConNet& cn;
	Fcon f;
	Variable overR;// gives 1/R
	std::vector<Variable> coords;
	Function * dist;
	std::vector<std::vector<Function*> > diff;
	int matsize;
	gsl_matrix * dynmat;
public:
	SphDynMatHarmOsc(Torus<Sphere>& bbox,HarmPot<Torus<Sphere> >& ppot,
			SphereConNet& ccn);
	~SphDynMatHarmOsc(){gsl_matrix_free(dynmat);}
	void calc_mat();//finds the dynamical matrix
	gsl_matrix * get_mat(){return dynmat;}//retrieves dynamical matrix
};

#endif /* SPHDYNMATHARMOSC_H_ */
