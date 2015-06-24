/*
 * SphereConNet.h
 *
 *  Created on: Jun 11, 2015
 *      Author: izzhov
 */

#ifndef SPHERECONNET_H_
#define SPHERECONNET_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Torus.h"
#include "Sphere.h"
#include "mm.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <cstddef>
#include <vector>

//this class is made for analysis at the *end*
class SphereConNet {
private:
	Torus<Sphere>& box;
	std::vector<std::vector<std::vector<int> > > con;
	std::vector<std::vector<gsl_vector*> > rs;
public:
	//bool is whether to auto-remove floaters
	SphereConNet(Torus<Sphere>& box, bool remove_em);

	double num_contacts();//number o' contacts
	//'Maxwell's criterion' for polydisperse spheres in PBC w/randomized packins
	int desired_contacts();
	int num_floaters();//number o' floaters
	//NOTE: floater functions only work if you've remove 'em

	//returns a vector listing which # particles are floaters
	//note indices start from zero
	std::vector<int> which_floaters();

	//removes all floaters (recursively)
	void rm_floaters();
	//removes one iteration of floaters, returns whether any floaters were removed
	bool rm_floaters_once();
	//precondition: part has any contacts
	bool is_floater(int i);
	bool is_full_3_floater(int i, int n); //precon: n>=2, for 2D
	bool is_full_4_floater(int i, int n); //precon: n>=3, for 3D
	bool is_3_floater(int x, int i, int j, int k);//x is for index, *2D* or *3D planar*
	bool is_4_floater(int x, int i, int j, int k, int l);//x=index, this is for *3D*
	//removes all the floater's contacts
	void rm(int i);

	//number of contacts on element i
	int num_contacts(int i);
	//returns con[i][j][k]
	int con_elem(int i, int j, int k);

	void print_con();//prints contacts to terminal
	void output_con(std::ofstream& dstream);//makes a file with all the pairs of points for the contacts

};

#endif /* SPHERECONNET_H_ */
