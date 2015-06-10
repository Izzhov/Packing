/*
 * HarmPot.h
 *
 *  Created on: Jun 9, 2015
 *      Author: izzhov
 */

#ifndef HARMPOT_H_
#define HARMPOT_H_

#include "Torus.h"
#include "mm.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <cmath>
#include <vector>

template <class B>//for Box
class HarmPot {
private:
	double del;
	double kk; double P;
	B& box;
	double U; std::vector<gsl_vector*> DU; //DU is for [0,1) coords
	std::vector<double> Dtaun; std::vector<gsl_vector*> Dtauv;
	double DL;
	std::vector<gsl_vector*> D_u; std::vector<gsl_vector*> D_v;
public:
	HarmPot(B& bbox);

	void set_P(double P);//re-calculates energy etc after setting
	void set_DU(int i, int j, double x);
	void set_Dtaun(int i, double x);
	void set_Dtauv(int i, int j, double x);
	void set_D_u(int i, int j, double x);
	void set_D_v(int i, int j, double x);

	double get_P(){return P;}
	double get_U(){return U;}
	gsl_vector * get_DU(int i){
		return DU[i];
	}
	gsl_vector * get_Dtauv(int i){
		return Dtauv[i];
	}
	double get_DU(int i, int j){
		return gsl_vector_get(DU[i],j);
	}
	double get_Dtaun(int i){
		return Dtaun[i];
	}
	double get_Dtauv(int i, int j){
		return gsl_vector_get(Dtauv[i],j);
	}
	double get_D_u(int i, int j){
		return gsl_vector_get(D_u[i],j);
	}
	double get_D_v(int i, int j){
		return gsl_vector_get(D_v[i],j);
	}
	double get_DL(){return DL;}

	//calculates energy and fills relevant variables
	//note these are *positive* derivatives
	void calc_U_F();
	void calc_U_F(int i, int j, int k);
	double dUdl(double lijk, double Rij);

	bool is_done(); //check if the total force is < d*avg_f
	//and total torque < d*d2 (10^(-2))*D(max dist of force appl)
	//and actual pressure is < twice & > half assigned press
	//see notebook (Press) and 6-9-15 for deets
	//assumes particle sizes are of similar order
	//assumes number of contacts is of order N

	double pressure();//calculates pressure
};

template<class B>
inline HarmPot<B>::HarmPot(B& bbox):box(bbox){
	del = std::pow(10,-6);//for is_done
	kk=1; P=std::pow(10,-4);//default values
	U=0; DL=0;
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextDU = gsl_vector_alloc(box.get_dim());
		gsl_vector * nextDtauv = gsl_vector_alloc(box.get_dim());
		gsl_vector * nextD_u = gsl_vector_alloc(box.get_dim());
		gsl_vector * nextD_v = gsl_vector_alloc(box.get_dim());
		gsl_vector_set_zero(nextDU); gsl_vector_set_zero(nextDtauv);
		gsl_vector_set_zero(nextD_u); gsl_vector_set_zero(nextD_v);
		DU.push_back(nextDU); Dtauv.push_back(nextDtauv);
		D_u.push_back(nextD_u); D_v.push_back(nextD_v);
		Dtaun.push_back(0);
	}
	calc_U_F();
}

template<class B>
inline void HarmPot<B>::set_P(double P){
	this->P = P;
	calc_U_F();
}

template<class B>
inline void HarmPot<B>::set_DU(int i, int j, double x) {
	gsl_vector_set(DU[i],j,x);
}

template<class B>
inline void HarmPot<B>::set_Dtaun(int i, double x) {
	Dtaun.at(i) = x;
}

template<class B>
inline void HarmPot<B>::set_Dtauv(int i, int j, double x) {
	gsl_vector_set(Dtauv[i],j,x);
}

template<class B>
inline void HarmPot<B>::set_D_u(int i, int j, double x){
	gsl_vector_set(D_u[i],j,x);
}

template<class B>
inline void HarmPot<B>::set_D_v(int i, int j, double x){
	gsl_vector_set(D_v[i],j,x);
}

template<class B>
inline void HarmPot<B>::calc_U_F() {
	double L = box.get_L(); int dim = box.get_dim();
	U = P*std::pow(L,dim); DL = dim*P*std::pow(L,dim-1);
	for(int i=0; i<box.get_N(); i++){
		set_Dtaun(i,0);
		for(int d=0; d<dim; d++){//initialize new forces to zero
			set_DU(i,d,0); set_Dtauv(i,d,0);
			set_D_u(i,d,0); set_D_v(i,d,0);
		}
		for(int j=0; j<box.get_N(); j++)
			for(int k=1; k<=mm::int_pow(2,dim);k++) calc_U_F(i,j,k);
	}
}

template<class B>
inline void HarmPot<B>::calc_U_F(int i, int j, int k) {
	if(i==j) return;
	double Rij = box.R(i,j); double lijk2 = box.ell2(i,j,k);
	if(lijk2>(Rij*Rij)) return;
	double lijk = sqrt(lijk2);
	if(i<j){//store energy
		U+=0.5*kk*(1-(lijk/Rij))*(1-(lijk/Rij));
	}
	double dudl = dUdl(lijk, Rij); double L = box.get_L();
	gsl_vector * f = box.ell_vec(i,j,k);//force direction (neg bc deriv)
	gsl_vector_scale(f,dudl/(L*lijk)); //force magnitude is good for [0,1)
	gsl_vector_add(DU[i],f);//adds this force to DU
	DL-=dudl/L;
	// TODO add torque rotations stuff
	gsl_vector_free(f);
}

template<class B>
inline double HarmPot<B>::dUdl(double lijk, double Rij) {
	return (kk/Rij)*(1-(lijk/Rij))*lijk;
}

template<class B>
inline bool HarmPot<B>::is_done() {
	bool id = true; int N = box.get_N(); int dim = box.get_dim();
	double L = box.get_L(); double fmr = (del*P*dim*std::pow(L,dim))/(double)N;
	for(int i=0; i<N; i++){
		double force; gsl_blas_ddot(DU[i],DU[i],&force);
		force = std::sqrt(force);
		if(force >= (fmr/(2*box.get_r(i)))){
			id = false; break;
		}
	}
	// TODO add torque check stuff
	return (id && -dim*std::pow(L,dim-1)*P<DL && DL<0.5*dim*std::pow(L,dim-1)*P);
}

template<class B>
inline double HarmPot<B>::pressure() {
	double L = box.get_L(); int dim = box.get_dim();
	return P-DL/(dim*std::pow(L,dim-1));
}

#endif /* HARMPOT_H_ */
