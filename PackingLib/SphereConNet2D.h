/*
 * SphereConNet2D.h
 *
 *  Created on: Jun 2, 2015
 *      Author: izzhov
 */

#ifndef SPHERECONNET2D_H_
#define SPHERECONNET2D_H_

#include "Torus2D.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <cstddef>
#include <vector>

//this class is made for analysis at the *end*
class SphereConNet2D {
private:
	Torus2D& torus;
	std::vector<std::vector<std::vector<int> > > con;
	std::vector<std::vector<gsl_vector*> > rs;
public:
	SphereConNet2D(Torus2D& torus);

	double num_contacts();//number o' contacts
	int desired_contacts();//2N-1, where N=# of non-floaters
	int num_floaters();//number o' floaters

	//returns a vector listing which # particles are floaters
	//note indices start from zero
	std::vector<int> which_floaters();

	//removes all floaters (recursively)
	void rm_floaters();
	//removes one iteration of floaters, returns whether any floaters were removed
	bool rm_floaters_once();
	//precondition: part has any contacts
	bool is_floater(int i);
	//removes all the floater's contacts
	void rm(int i);

	void print_con();
};

#endif /* SPHERECONNET2D_H_ */
