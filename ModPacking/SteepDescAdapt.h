/*
 * SteepDescAdapt.h
 *
 *  Created on: Jul 8, 2015
 *      Author: izzhov
 */

#ifndef STEEPDESCADAPT_H_
#define STEEPDESCADAPT_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include "mm.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <vector>
#include <cmath>

template <class P, class B>//P is potential, B is box
class SteepDescAdapt {
private:
	double eps; double epsL;
	std::vector<double> press;
	P& pot;
	B& box;
public:
	SteepDescAdapt(P& ppot, B& bbox);
	void minimize();
	void minimize(int n);//minimizes at successively smaller pressures
	bool exp_minimize(int n);//minimizes with P=10^(-4)*0.9^i, i<n
	void next();
	void move();
	void translate();
	void rotate();

	void calc_eps();//P/F*0.01 (note: assumes pressure is 10^-4 or less)
	double max_force();//max of norm of each DU * L
};

template<class P, class B>
inline SteepDescAdapt<P, B>::SteepDescAdapt(P& ppot, B& bbox):pot(ppot),box(bbox){
	eps=0.1; epsL=0.1;
	for(int b=-4; b>-20; b=b-2){
		press.push_back(std::pow(10,b));
	}
}

template<class P, class B>
inline void SteepDescAdapt<P, B>::minimize() {
	int i=0;
	while(true){
		next(); i++;
		if(i%1000==0){
			if(pot.is_done()) break;
			if(i>=1000000) break;
			calc_eps();
		}
	}
	if(i>2000){
		std::cout << i << " " << "flt: " << pot.is_float() << std::endl;
	}
}

template<class P, class B>
inline void SteepDescAdapt<P, B>::minimize(int n) {
	for(int i=0; i<n; i++){
		pot.set_P(press[i]);
		minimize();
	}
}

template<class P, class B>
inline bool SteepDescAdapt<P, B>::exp_minimize(int n){
	for(int i=0; i<n; i++){
		minimize();
		pot.set_P(pot.get_P()*0.9);
	}
	return pot.is_float();
}

template<class P, class B>
inline void SteepDescAdapt<P, B>::next() {
	move(); pot.calc_U_F();
}

template<class P, class B>
inline void SteepDescAdapt<P, B>::move() {
	translate();
	if(box.sym()>1) rotate();
}

template<class P, class B>
inline void SteepDescAdapt<P, B>::translate() {
	double L = box.get_L(); double DL = pot.get_DL();
	double prs = pot.get_P();
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_scale(pot.get_DU(i),eps);
		//now cap each component at P
		for(int d=0; d<box.get_dim(); d++){
			if(std::abs(gsl_vector_get(pot.get_DU(i),d))>0.01*prs)
				gsl_vector_set(pot.get_DU(i),d,
						mm::sign(0.01*prs,gsl_vector_get(pot.get_DU(i),d)));
		}
		gsl_vector_sub(box.get_1_pos(i),pot.get_DU(i));
		box.modify(i);//makes it all mod 1
	}
	box.set_L(L-epsL*DL);
}

template<class P, class B>
inline void SteepDescAdapt<P, B>::rotate() {
	double prs = pot.get_P();
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_scale(pot.get_D_u(i),eps);
		//now cap each component at P
		for(int d=0; d<box.get_dim(); d++){
			if(std::abs(gsl_vector_get(pot.get_D_u(i),d))>0.01*prs)
				gsl_vector_set(pot.get_D_u(i),d,
						mm::sign(0.01*prs,gsl_vector_get(pot.get_D_u(i),d)));
		}
		gsl_vector_sub(box.get_u(i),pot.get_D_u(i));
		if(box.sym()==3){
			gsl_vector_scale(pot.get_D_v(i),eps);
			//now cap each component at P
			for(int d=0; d<box.get_dim(); d++){
				if(std::abs(gsl_vector_get(pot.get_D_v(i),d))>0.01*prs)
					gsl_vector_set(pot.get_D_v(i),d,
							mm::sign(0.01*prs,gsl_vector_get(pot.get_D_v(i),d)));
			}
			gsl_vector_sub(box.get_v(i),pot.get_D_v(i));
		}
		box.normalize(i);
	}
}

template<class P, class B>
inline void SteepDescAdapt<P, B>::calc_eps() {
	double mf = max_force(); double prs = pot.get_P();
	if(mf<std::pow(10,-13)) eps = 0.01*prs/std::pow(10.0,-12.0);
	else{
		double newep=0.01*prs/mf;
		if(newep>0.1){
			if(newep<std::pow(10,6)) eps=newep;
			else eps=std::pow(10,6);
		}
		else eps=0.1;
	}
}

template<class P, class B>
inline double SteepDescAdapt<P, B>::max_force() {
	double mf = 0; double cf = 0;
	for(int i=0; i<box.get_N(); i++)
		if((cf=box.get_L()*gsl_blas_dnrm2(pot.get_DU(i)))>mf) mf=cf;
	return mf;
}

#endif /* STEEPDESCADAPT_H_ */
