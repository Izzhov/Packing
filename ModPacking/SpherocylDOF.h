/*
 * SpherocylDOF.h
 *
 *  Created on: Jul 6, 2015
 *      Author: izzhov
 */

#ifndef SPHEROCYLDOF_H_
#define SPHEROCYLDOF_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include "mm.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <cmath>
#include <vector>

class SpherocylDOF {
private:
	gsl_vector * pos;
	gsl_vector * u;
	double r;
	double a;
public:
	SpherocylDOF();
	SpherocylDOF(gsl_vector * pos, std::vector<gsl_vector*> size);//blank
	SpherocylDOF(gsl_vector * pos, gsl_vector * u, std::vector<gsl_vector*> size);
	SpherocylDOF(gsl_vector * pos, gsl_vector * u, gsl_vector * v, std::vector<gsl_vector*> size);//blank
	//makes shell copy that's ratio times as big (see 7-9-15 pg 2)
	SpherocylDOF(SpherocylDOF * orig, double ratio);
	/*
	~SpherocylDOF();
	//copy constructor so vector copies copy content not pointers
	SpherocylDOF(const SpherocylDOF& other);
	*/

	static int sym(int dim){return 2;}//always has 1 rot. vec
	static int ncon(int dim){return dim;}//# of possible pts of contact

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
	double lsl();//longest shape length

	gsl_vector * lisljs(SpherocylDOF s, int k, double L, int ncon);//finds lis and ljs from algorithm
	gsl_vector * F_loc(SpherocylDOF s, int k, double L, int ncon);//place at which force is applied
	gsl_vector * ell_vec(SpherocylDOF s, int k, double L, int ncon);//in [0,1) space
	//also in [0,1) space
	double ell2(SpherocylDOF s, int k, double L, int ncon);//dist btwn in k-th quadrant (see 6-8-15)
	//now, versions of all these which are quadrantless:
	gsl_vector * lisljs(SpherocylDOF s, double L, int ncon);
	gsl_vector * F_loc(SpherocylDOF s, double L, int ncon);
	gsl_vector * ell_vec(SpherocylDOF s, double L, int ncon);
	double ell2(SpherocylDOF s, double L, int ncon);

	bool touch(SpherocylDOF s, double L);//checks if they touch
	bool is_far(SpherocylDOF s, double L);//checks if outside nbrlist-zone

	//sets u to magnitude 1
	void normalize();

	//does nothing for spherocyls
	gsl_vector * get_v();
};

#endif /* SPHEROCYLDOF_H_ */
