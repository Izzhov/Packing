/*
 * Sphere.h
 *
 *  Created on: Jun 8, 2015
 *      Author: izzhov
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include "mm.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <cmath>
#include <vector>

class Sphere {
private:
	gsl_vector * pos;
	double r;
public:
	Sphere();
	Sphere(gsl_vector * pos, std::vector<gsl_vector*> size);
	//these 2 are blank since spheres don't need u and v
	Sphere(gsl_vector * pos, gsl_vector * u, std::vector<gsl_vector*> size);
	Sphere(gsl_vector * pos, gsl_vector * u, gsl_vector * v, std::vector<gsl_vector*> size);

	static int sym(int dim){return 1;}//sphere doesn't need rotation vecs

	void set_pos(gsl_vector * pos);
	void set_r(double r);
	void set_pos_coord(int i, double x);//sets i-th coord to x; mods 1 for torus

	gsl_vector * get_pos();
	double get_r();
	double get_pos_coord(int i);

	double volume();//volume or area depending on dimension

	gsl_vector * ell_vec(Sphere s, int k, double L);
	double ell2(Sphere s, int k, double L);//dist btwn in k-th quadrant (see 6-8-15)

	//does nothing for spheres
	void normalize();
	double max_d();
	gsl_vector * get_u();
	gsl_vector * get_v();
	double I();
	gsl_vector * F_loc(Sphere s, int k, double L);
};

#endif /* SPHERE_H_ */
