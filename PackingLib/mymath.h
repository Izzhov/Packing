/*
 * mymath.h
 *
 *  Created on: May 22, 2015
 *      Author: root
 */

#ifndef MYMATH_H_
#define MYMATH_H_

#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

class mymath {
public:
	static double pmod(double x, double y){
		return fmod(y+fmod(x,y),y);
	}

	static double dist2d(double x1, double y1, double x2, double y2){
		return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	}

	static double sign(double a, double b){
		if (b>=0) return std::abs(a);
		else return -std::abs(a);
	}
	//stores the cross product of u and v (3-vectors) in product
	static void cross_product(const gsl_vector * u, const gsl_vector * v, gsl_vector * product){
		double p1 = gsl_vector_get(u,1)*gsl_vector_get(v,2)-
				gsl_vector_get(u,2)*gsl_vector_get(v,1);
		double p2 = gsl_vector_get(u,2)*gsl_vector_get(v,0)-
				gsl_vector_get(u,0)*gsl_vector_get(v,2);
		double p3 = gsl_vector_get(u,0)*gsl_vector_get(v,1)-
				gsl_vector_get(u,1)*gsl_vector_get(v,0);
		gsl_vector_set(product,0,p1);
		gsl_vector_set(product,1,p2);
		gsl_vector_set(product,2,p3);
	}

	static int sgn(double x){
		if(x>0) return 1;
		else if(x==0) return 0;
		else return -1;
	}
	//takes cross product of (z-hat) X (2d planar vector)
	static void cross2d(const gsl_vector * u, gsl_vector * product){
		double p1 = -gsl_vector_get(u,1);
		double p2 = gsl_vector_get(u,0);
		gsl_vector_set(product,0,p1);
		gsl_vector_set(product,1,p2);
	}

};

#endif /* MYMATH_H_ */
