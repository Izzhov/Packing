/*
 * FullSecSteepDesc.h
 *
 *  Created on: Jul 13, 2015
 *      Author: izzhov
 */

#ifndef FULLSECSTEEPDESC_H_
#define FULLSECSTEEPDESC_H_

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
class FullSecSteepDesc {
private:
	double del; double del0;
	double alph; double eps; double seceps;
	double sig; int imax; int jmax; //see conjgrad paper for meaning of these
	P& pot; B& box;
	double expmin;//fraction do decrease pressure and eps by each exp_min
public:
	FullSecSteepDesc(P& ppot, B& bbox);
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

	//how accurate the simulation really got (sqrt(del/del0))
	double eps_reached(){return std::sqrt(del/del0);}
};

template<class P, class B>
inline FullSecSteepDesc<P, B>::FullSecSteepDesc(P& ppot, B& bbox):pot(ppot),box(bbox){
	alph=0; eps=std::pow(10,-4); seceps=std::pow(10,-4);
	sig=0.1; imax = 1000000; jmax = 1000;
	//std::cout << "Seceps: " << seceps << std::endl;
	expmin = 0.5;
}

template<class P, class B>
inline bool FullSecSteepDesc<P, B>::minimize() {
	int i=0;
	double eta; double etaprev;
	//r stores negative deriv; rsec stores pos deriv after moving sig0
	//rold stores the descent direction always used in the do-while loop
	gsl_vector * r = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	gsl_vector * rold = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	gsl_vector * rsec = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	//r=b-Ax
	pot.calc_U_F();
	store_deriv(r,-1);
	//del=rTr
	del = unit_fix(r); //formerly gsl_blas_ddot(r,r,&del);
	del0=del;
	while(i<imax && del>eps*eps*del0){
		int j=0; store_deriv(rold,-1);
		alph=-sig;
		move(sig); pot.calc_U_F(); store_deriv(rsec,1);//find f'(x+sig*d)
		gsl_blas_ddot(rsec,r,&etaprev);
		//old location of the move(-sig)... stuff:
		move(-sig); load_deriv(r,-1); //same desc-dir
		double alphbd = 100;//max 'n' min for desc dir
		//Uprev: double Uprev = pot.get_U();
		//have max and min been hit for current alphbd?
		bool amax = false; bool amin = false;
		do{
			gsl_blas_ddot(r,rold,&eta); eta=-eta;
			if(mm::sgn(eta)!=mm::sgn(alph*eta/(etaprev-eta)))
				alph = alph*eta/(etaprev-eta);
			else alph=-2*mm::sgn(eta)*std::abs(alph);
			/*
			if(j==0){
				//if(alph<0) break;
				move(-sig); load_deriv(r,-1); //same desc-dir
			}
			*/
			if(amax && amin){alphbd=alphbd/2.0; amax=false; amin=false;}
			if(alph>alphbd){
				alph=alphbd;
				amax=true;
			}
			else if(alph<-alphbd){
				alph=-alphbd;
				amin=true;
			}
			else{alphbd=std::abs(alph);}
			load_deriv(rold,-1); //deploy desc-dir before moving
			move(alph); pot.calc_U_F(); store_deriv(r,-1);
			etaprev=eta;
			/*Uprev:
			if(pot.get_U()<Uprev) Uprev = pot.get_U();
			else{
				load_deriv(rold,-1); //deploy desc-dir before moving
				move(-alph); pot.calc_U_F(); store_deriv(r,-1); break;
			}
			*/
			j++;
		}while(j<jmax && alph*alph*del>seceps*seceps);
		/*
		if(i%10000==0 || j>100){
			std::cout << "i: " << i << " j: " << j << " U: " << pot.get_U() << std::endl;
		}
		*/
		del = unit_fix(r);
		del0 = box.get_dim()*box.get_dim()*pot.get_P()*pot.get_P()*std::pow(box.get_L(),2*box.get_dim()-2);
		/*
		if(i%10000==0){
			std::cout << i << std::endl;
			//std::cout << "L: " << box.get_L() << std::endl;
			//std::cout << "L^2? (int): " << std::pow(box.get_L(),2*box.get_dim()-2) << std::endl;
			//std::cout << "(2PL)^2: " << (2*pot.get_P()*box.get_L())*(2*pot.get_P()*box.get_L()) << std::endl;
			std::cout << "del0: " << del0 << std::endl;
			std::cout << "del: " << del << std::endl;
		}
		*/
		i++;
	}
	gsl_vector_free(r); gsl_vector_free(rsec);
	std::cout << i << " out of " << imax << std::endl;
	return i<imax;
}

template<class P, class B>
inline bool FullSecSteepDesc<P, B>::exp_minimize(int n) {
	bool madeit = false;
	for(int i=0; i<n; i++){
		madeit = minimize();
		pot.set_P(pot.get_P()*expmin);
		eps*=expmin;
	}
	return madeit;
}

template<class P, class B>
inline void FullSecSteepDesc<P, B>::move(double alphsig) {
	translate(alphsig);
	if(box.sym()>1) rotate(alphsig);
}

template<class P, class B>
inline void FullSecSteepDesc<P, B>::translate(double alphsig) {
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
	//double Lexp = -alphsig*L*DL;
	box.set_L(L-alphsig*DL);
}

template<class P, class B>
inline void FullSecSteepDesc<P, B>::rotate(double alphsig) {
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
inline void FullSecSteepDesc<P, B>::store_deriv(gsl_vector* r, int scale) {
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
inline void FullSecSteepDesc<P, B>::load_deriv(gsl_vector* r, int scale) {
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
inline double FullSecSteepDesc<P, B>::unit_fix(gsl_vector* r) {
	double dotprod=0; double component;
	for(int i=0; i<r->size; i++){
		component = gsl_vector_get(r,i);
		if(i%box.sym()==0 && i<(r->size)-1) component*=box.get_L();
		dotprod+=component*component;
	}
	return dotprod;
}

#endif /* FULLSECSTEEPDESC_H_ */
