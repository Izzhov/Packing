/*
 * SecSteepDesc.h
 *
 *  Created on: Jul 8, 2015
 *      Author: izzhov
 */

#ifndef SECSTEEPDESC_H_
#define SECSTEEPDESC_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <sstream>

#include "mm.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <vector>
#include <cmath>

template <class P, class B>//P is potential, B is box
class SecSteepDesc {
private:
	double alph; double eps;
	double sig; int imax; //see conjgrad paper for meaning of these
	P& pot; B& box;
public:
	SecSteepDesc(P& ppot, B& bbox);
	bool minimize();
	bool exp_minimize(int n);//minimizes with P=10^(-4)*0.9^i, i<n

	void move(double alphsig);
	void translate(double alphsig);
	void rotate(double alphsig);

	void store_deriv(gsl_vector * r, int scale);//scale to store pos or neg
	void load_deriv(gsl_vector * r, int scale);
	//returns the vector^2 but with every sym-th element mult'd by L
	//(to fix the units on the forces - see7-9-15 whiteback)
	double unit_fix(gsl_vector * r);
};

template<class P, class B>
inline SecSteepDesc<P, B>::SecSteepDesc(P& ppot, B& bbox):pot(ppot),box(bbox){
	alph=0; eps=std::pow(10,-4); sig=0.1; imax = 1000000;
}

template<class P, class B>
inline bool SecSteepDesc<P, B>::minimize() {
	int i=0; double del; double del0;
	double eta; double etaprev;
	//r stores negative deriv; rsec stores pos deriv after moving sig0
	gsl_vector * r = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	gsl_vector * rsec = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	//r=b-Ax
	pot.calc_U_F();
	store_deriv(r,-1);
	//del=rTr
	del = unit_fix(r); //formerly gsl_blas_ddot(r,r,&del);
	del0=del;
	while(i<imax && del>eps*eps*del0){
		alph=-sig;
		move(sig); pot.calc_U_F(); store_deriv(rsec,1);//find f'(x+sig*d)
		gsl_blas_ddot(rsec,r,&etaprev);
		gsl_blas_ddot(r,r,&eta); eta=-eta;
		alph = alph*eta/(etaprev-eta);
		if(alph>100) alph=100;
		if(alph>=sig){//only do second move if it's in the right direction
			//need to load old derivs before second move
			load_deriv(r,-1);
			move(alph-sig); pot.calc_U_F(); store_deriv(r,-1);
		}
		del = unit_fix(r);
		del0 = box.get_dim()*box.get_dim()*pot.get_P()*pot.get_P()*std::pow(box.get_L(),2*box.get_dim()-2);
		if(i%10000==0){
			std::cout << i << std::endl;
			//std::cout << "L: " << box.get_L() << std::endl;
			//std::cout << "L^2? (int): " << std::pow(box.get_L(),2*box.get_dim()-2) << std::endl;
			//std::cout << "(2PL)^2: " << (2*pot.get_P()*box.get_L())*(2*pot.get_P()*box.get_L()) << std::endl;
			std::cout << "del0: " << del0 << std::endl;
			std::cout << "del: " << del << std::endl;
		}
		i++;
	}
	gsl_vector_free(r); gsl_vector_free(rsec);
	std::cout << i << " out of " << imax << std::endl;
	return i>=imax;
}

template<class P, class B>
inline bool SecSteepDesc<P, B>::exp_minimize(int n) {
	bool madeit = false;
	for(int i=0; i<n; i++){
		madeit = minimize();
		pot.set_P(pot.get_P()*0.9);
		eps*=0.9;
	}
	return madeit;
}

template<class P, class B>
inline void SecSteepDesc<P, B>::move(double alphsig) {
	translate(alphsig);
	if(box.sym()>1) rotate(alphsig);
}

template<class P, class B>
inline void SecSteepDesc<P, B>::translate(double alphsig) {
	double L = box.get_L(); double DL = pot.get_DL();
	//double moveL;//to figure out how to move L (see 7-8-15 page 1);
	gsl_vector * mover = gsl_vector_alloc(box.get_dim());
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_memcpy(mover,pot.get_DU(i));
		gsl_vector_scale(mover,alphsig);
		gsl_vector_sub(box.get_1_pos(i),mover);
		box.modify(i);//makes it all mod 1
	}
	gsl_vector_free(mover);
	//if(DL>=0) moveL=std::min(alphsig*DL,0.1*L);
	//else moveL=std::max(alphsig*DL,alphL*DL);
	box.set_L(L-alphsig*DL);
}

template<class P, class B>
inline void SecSteepDesc<P, B>::rotate(double alphsig) {
	gsl_vector * mover = gsl_vector_alloc(box.get_dim());
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_memcpy(mover,pot.get_D_u(i));
		gsl_vector_scale(mover,alphsig);
		gsl_vector_sub(box.get_u(i),mover);
		if(box.sym()==3){
			gsl_vector_memcpy(mover,pot.get_D_v(i));
			gsl_vector_scale(mover,alphsig);
			gsl_vector_sub(box.get_v(i),mover);
		}
		box.normalize(i);
	}
	gsl_vector_free(mover);
}

template<class P, class B>
inline void SecSteepDesc<P, B>::store_deriv(gsl_vector* r, int scale) {
	int i=0;
	for(int j=0; j<box.get_N(); j++){
		for(int d=0; d<box.get_dim(); d++){
			gsl_vector_set(r,i,scale*gsl_vector_get(pot.get_DU(j),d));
			i++;
			if(box.sym()>1){
				gsl_vector_set(r,i,scale*gsl_vector_get(pot.get_D_u(j),d));
				i++;
				if(box.sym()==3){
					gsl_vector_set(r,i,scale*gsl_vector_get(pot.get_D_v(j),d));
					i++;
				}
			}
		}
	}
	gsl_vector_set(r,i,scale*pot.get_DL());
}

template<class P, class B>
inline void SecSteepDesc<P, B>::load_deriv(gsl_vector* r, int scale) {
	int i=0;
	for(int j=0; j<box.get_N(); j++){
		for(int d=0; d<box.get_dim(); d++){
			pot.set_DU(j,d,scale*gsl_vector_get(r,i));
			i++;
			if(box.sym()>1){
				pot.set_D_u(j,d,scale*gsl_vector_get(r,i));
				i++;
				if(box.sym()==3){
					pot.set_D_v(j,d,scale*gsl_vector_get(r,i));
					i++;
				}
			}
		}
	}
	pot.set_DL(scale*gsl_vector_get(r,i));
}

template<class P, class B>
inline double SecSteepDesc<P, B>::unit_fix(gsl_vector* r) {
	double dotprod=0; double component;
	for(int i=0; i<r->size; i++){
		component = gsl_vector_get(r,i);
		if(i%box.sym()==0 && i<(r->size)-1) component*=box.get_L();
		dotprod+=component*component;
	}
	return dotprod;
}

#endif /* SECSTEEPDESC_H_ */
