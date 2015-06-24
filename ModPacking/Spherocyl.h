/*
 * Spherocyl.h
 *
 *  Created on: Jun 16, 2015
 *      Author: izzhov
 */

#ifndef SPHEROCYL_H_
#define SPHEROCYL_H_

#include "mm.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <cmath>
#include <vector>

class Spherocyl {
private:
	gsl_vector * pos;
	gsl_vector * u;
	double r;
	double a;
public:
	Spherocyl();
	Spherocyl(gsl_vector * pos, std::vector<gsl_vector*> size);//blank
	Spherocyl(gsl_vector * pos, gsl_vector * u, std::vector<gsl_vector*> size);
	Spherocyl(gsl_vector * pos, gsl_vector * u, gsl_vector * v, std::vector<gsl_vector*> size);//blank

	static int sym(int dim){return 2;}//always has 1 rot. vec

	void set_pos(gsl_vector * pos);
	void set_u(gsl_vector * u);
	void set_r(double r);
	void set_a(double a);
	void set_pos_coord(int i, double x);//sets i-th coord to x; mods 1 for torus
	void set_u_coord(int i, double x);

	gsl_vector * get_pos(){return pos;}
	gsl_vector * get_u(){return u;}
	double get_r(){return r;}
	double get_a(){return a;}
	double get_pos_coord(int i){return gsl_vector_get(pos,i);}
	double get_u_coord(int i){return gsl_vector_get(u,i);}

	double sphere_volume();
	double cyl_volume();
	double volume();//volume or area depending on dimension
	double I();//moment of inertia
	double max_d();//largest distance from center that force can be applied

	gsl_vector * lisljs(Spherocyl s, int k, double L);//finds lis and ljs from algorithm
	gsl_vector * F_loc(Spherocyl s, int k, double L);//place at which force is applied
	gsl_vector * ell_vec(Spherocyl s, int k, double L);
	double ell2(Spherocyl s, int k, double L);//dist btwn in k-th quadrant (see 6-8-15)

	//sets u to magnitude 1
	void normalize();
};

#endif /* SPHEROCYL_H_ */
