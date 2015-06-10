/*
 * mm.h
 *
 *  Created on: Jun 8, 2015
 *      Author: izzhov
 */

#ifndef MM_H_
#define MM_H_

#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

class mm {
public:
	static double pmod(double x, double y){
		return fmod(y+fmod(x,y),y);
	}
	static int int_pow(int x, int p){//NOTE: ONLY FOR NON-NEGATIVE INTS
		if (p==0) return 1;
		if (p==1) return x;

		int tmp = int_pow(x,p/2);
		if (p%2==0) return tmp*tmp;
		else return x*tmp*tmp;
	}
	//make sure it's mod 1
	static gsl_vector* rel(const gsl_vector* v1, const gsl_vector* v2, int k){
		gsl_vector* result = gsl_vector_alloc(v1->size);
		double coord; //coordinate

		for (unsigned int i=0; i<(v1->size); i++){
			coord = pmod(gsl_vector_get(v2,i)-gsl_vector_get(v1,i),1);
			int minq = int_pow(2,i);//minimum value of (k mod 2*minq) to -=1
			if((k%(2*minq))>=minq) coord-=1;
			gsl_vector_set(result,i,coord);
		}
		return result;
	}
};

#endif /* MM_H_ */
