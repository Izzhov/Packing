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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <stdarg.h>

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
	//checks if 2 vectors are antiparallel
	static bool antiparallel(const gsl_vector* v1, const gsl_vector* v2){
		double normprod = gsl_blas_dnrm2(v1)*gsl_blas_dnrm2(v2);
		double dotprod; gsl_blas_ddot(v1,v2,&dotprod);
		double prodsum = normprod+dotprod;//if small->antiparallel
		return (prodsum<std::pow(10,-14));
	}
	static int sgn(double x){
		if(x>0) return 1;
		else if(x==0) return 0;
		else return -1;
	}
	static double det_get(gsl_matrix * A, int inPlace){
		/*
		 * inPlace == 1 => A is replaced with the LU decomposed copy.
		 * inPlace == 0 => A is retained, and a copy is used for LU.
		 */

		double det;
		int signum;
		gsl_permutation * p = gsl_permutation_alloc(A->size1);
		gsl_matrix * tmpA;

		if (inPlace)
			tmpA=A;
		else{
			tmpA = gsl_matrix_alloc(A->size1, A->size2);
			gsl_matrix_memcpy(tmpA,A);
		}

		gsl_linalg_LU_decomp(tmpA,p,&signum);
		det = gsl_linalg_LU_det(tmpA,signum);
		gsl_permutation_free(p);
		if(!inPlace) gsl_matrix_free(tmpA);

		return det;
	}
	static gsl_vector* cross(int n_args, ...){//vectors all have length n_args+1
		gsl_vector * crossprod = gsl_vector_alloc(n_args+1);
		gsl_matrix * A = gsl_matrix_calloc(n_args+1,n_args+1);
		gsl_vector * topaste;//keep track of what you're pasting
		//now paste the vectors into the first n_args rows of A
		va_list ap;
		va_start(ap, n_args);
		for(int i=0; i<n_args; i++){
			topaste = va_arg(ap,gsl_vector*);
			gsl_matrix_set_row(A,i,topaste);
		}
		va_end(ap);
		//set determinants as cross product components (see 6-11-15 page 2)
		for(int i=0; i<(n_args+1); i++){
			gsl_matrix_set(A,n_args,i,1);
			gsl_vector_set(crossprod,i,det_get(A,0));
			gsl_matrix_set(A,n_args,i,0);
		}
		gsl_matrix_free(A);
		return crossprod;
	}
};

#endif /* MM_H_ */
