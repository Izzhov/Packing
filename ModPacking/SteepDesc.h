/*
 * SteepDesc.h
 *
 *  Created on: Jun 9, 2015
 *      Author: izzhov
 */

#ifndef STEEPDESC_H_
#define STEEPDESC_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <vector>
#include <cmath>

template <class P, class B>//P is potential, B is box
class SteepDesc {
private:
	double eps;
	std::vector<double> press;
	P& pot;
	B& box;
public:
	SteepDesc(P& ppot, B& bbox);
	void minimize();
	void minimize(int n);//minimizes at successively smaller pressures
	bool exp_minimize(int n);//minimizes with P=10^(-4)*0.9^i, i<n
	void next();
	void move();
	void translate();
	void rotate();
};

template<class P, class B>
inline SteepDesc<P, B>::SteepDesc(P& ppot, B& bbox):pot(ppot),box(bbox){
	eps=0.1;
	for(int b=-4; b>-20; b=b-2){
		press.push_back(std::pow(10,b));
	}
}

template<class P, class B>
inline void SteepDesc<P, B>::minimize() {
	int i=0;
	while(true){
		next(); i++;
		if(i%1000==0){
			if(pot.is_done()) break;
			if(i>=1000000) break;
		}
	}
	if(i>2000){
		std::cout << i << " " << "flt: " << pot.is_float() << std::endl;
	}
}

template<class P, class B>
inline void SteepDesc<P, B>::minimize(int n) {
	for(int i=0; i<n; i++){
		pot.set_P(press[i]);
		minimize();
	}
}

template<class P, class B>
inline bool SteepDesc<P, B>::exp_minimize(int n){
	for(int i=0; i<n; i++){
		minimize();
		pot.set_P(pot.get_P()*0.9);
	}
	return pot.is_float();
}

template<class P, class B>
inline void SteepDesc<P, B>::next() {
	move(); pot.calc_U_F();
}

template<class P, class B>
inline void SteepDesc<P, B>::move() {
	translate();
	if(box.sym()>1) rotate();
}

template<class P, class B>
inline void SteepDesc<P, B>::translate() {
	double L = box.get_L(); double DL = pot.get_DL();
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_scale(pot.get_DU(i),eps);
		gsl_vector_sub(box.get_1_pos(i),pot.get_DU(i));
		box.modify(i);//makes it all mod 1
	}
	box.set_L(L-eps*DL);
}

template<class P, class B>
inline void SteepDesc<P, B>::rotate() {
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_scale(pot.get_D_u(i),eps);
		gsl_vector_sub(box.get_u(i),pot.get_D_u(i));
		if(box.sym()==3){
			gsl_vector_scale(pot.get_D_v(i),eps);
			gsl_vector_sub(box.get_v(i),pot.get_D_v(i));
		}
		box.normalize(i);
	}
}

#endif /* STEEPDESC_H_ */
