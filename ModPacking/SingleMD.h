/*
 * SingleMD.h
 *
 *  Created on: Jun 22, 2015
 *      Author: izzhov
 */

#ifndef SINGLEMD_H_
#define SINGLEMD_H_

#include "mm.h"
#include "HarmPot.h"
#include "Torus.h"
#include "TimeRNGNormal.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <vector>
#include <cmath>

template <class P, class B>//P is potential, B is box
class SingleMD {
private:
	TimeRNGNormal& trng;
	double PE0; //minimal PE
	double avgU; double avgU2;//for finding avg and std-dev
	double kT; //10^(-8)
	B box; P * pot;
	std::vector<gsl_vector*> oldDU;
	std::vector<double> oldDtaun; std::vector<gsl_vector*> oldDtauv;
	std::vector<gsl_vector*> vs;
	std::vector<double> wns; std::vector<gsl_vector*> wvs;
	double delt; int numstep;
public:
	SingleMD(B bbox, double ddelt, TimeRNGNormal& ttrng);
	//copy constructor, ndef is so not mixed up with default copy
	SingleMD(SingleMD<P,B>& copy, int ndef);

	double get_PE0(){return PE0;}
	double get_avgU(){return avgU;}
	double get_avgU2(){return avgU2;}
	double get_kT(){return kT;}
	B get_box(){return box;}
	gsl_vector* get_oldDU(int i){return oldDU[i];}
	std::vector<double> get_oldDtaun(){return oldDtaun();}
	gsl_vector* get_oldDtauv(int i){return oldDtauv[i];}
	gsl_vector* get_vs(int i){return vs[i];}
	std::vector<double> get_wns(){return wns;}
	gsl_vector* get_wvs(int i){return wvs[i];}
	double get_delt(){return delt;}
	int get_numstep(){return numstep;}

	void populate_v();
	void populate_wn();
	void populate_wv();

	double energy();

	void next();

	void move_x();
	void rot();
	void store_DU();
	void store_Dtaun();
	void store_Dtauv();
	void move_v();
	void move_w();
};

template<class P, class B>
inline SingleMD<P, B>::SingleMD(B bbox, double ddelt, TimeRNGNormal& ttrng):trng(ttrng){
	box = bbox; delt = ddelt; numstep=1;
	pot = new P(box);
	PE0 = (*pot).get_U(); kT = std::pow(10,-8);
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextDU = gsl_vector_alloc(box.get_dim());
		gsl_vector * nextDtauv = gsl_vector_alloc(box.get_dim());
		gsl_vector_set_zero(nextDU); gsl_vector_set_zero(nextDtauv);
		oldDU.push_back(nextDU); oldDtauv.push_back(nextDtauv);
		oldDtaun.push_back(0);
	}
	populate_v();
	if(box.sym()>1){
		if(box.get_dim()==2) populate_wn();
		else populate_wv();
	}
	double nrg = energy();
	avgU = nrg; avgU2 = nrg*nrg;
}

template<class P, class B>
inline SingleMD<P, B>::SingleMD(SingleMD<P,B>& copy, int ndef){
	PE0 = copy.get_PE0(); avgU = copy.get_avgU(); avgU2 = copy.get_avgU2();
	kT = copy.get_kT(); oldDtaun = copy.get_oldDtaun(); wns = copy.get_wns();
	delt = copy.get_delt(); numstep = copy.get_numstep();
	B bbox(copy.get_box(),0); box = bbox;
	pot = new P(box);
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextDU = gsl_vector_alloc(box.get_dim());
		gsl_vector * nextDtauv = gsl_vector_alloc(box.get_dim());
		gsl_vector * nextv = gsl_vector_alloc(box.get_dim());
		gsl_vector_memcpy(nextDU,copy.get_oldDU(i));
		gsl_vector_memcpy(nextDtauv,copy.get_oldDtauv(i));
		gsl_vector_memcpy(nextv,copy.get_vs(i));
		oldDU.push_back(nextDU); oldDtauv.push_back(nextDtauv);
		vs.push_back(nextv);
	}
	if(box.sym()>1 && box.get_dim()==3) for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextwv = gsl_vector_alloc(box.get_dim());
		gsl_vector_memcpy(nextwv,copy.get_wvs(i));
		wvs.push_back(nextwv);
	}
}

template<class P, class B>
inline void SingleMD<P, B>::populate_v(){
	double skT = std::sqrt(kT);
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextv = gsl_vector_alloc(box.get_dim());
		for(int d=0; d<box.get_dim(); d++)
			gsl_vector_set(nextv,d,skT*trng.num());
		vs.push_back(nextv);
	}
}

template<class P, class B>
inline void SingleMD<P, B>::populate_wn(){
	double skT = std::sqrt(kT);
	for(int i=0; i<box.get_N(); i++){
		double sI = std::sqrt(box.I(i));
		wns.push_back((skT/sI)*trng.num());
	}
}

template<class P, class B>
inline void SingleMD<P, B>::populate_wv() {
	double skT = std::sqrt(kT);
	if(box.sym()==2) for(int i=0; i<box.get_N(); i++){
		double sI = std::sqrt(box.I(i));
		gsl_vector * nextwv = gsl_vector_calloc(box.get_dim());
		//need to make vectors for the 2 degrees of freedom
		gsl_vector * wv1 = mm::perp_3(box.get_u(i),1);
		gsl_vector * wv2 = mm::perp_3(box.get_u(i),2);
		mm::normalize(wv1); mm::normalize(wv2);
		gsl_vector_scale(wv1,(skT/sI)*trng.num());
		gsl_vector_scale(wv2,(skT/sI)*trng.num());
		gsl_vector_add(nextwv,wv1);
		gsl_vector_add(nextwv,wv2);
		wvs.push_back(nextwv);
		gsl_vector_free(wv1); gsl_vector_free(wv2);
	}
	else{
		//TODO work out populatin this with Imat
	}
}

template<class P, class B>
inline double SingleMD<P, B>::energy() {
	double nrg = (*pot).get_U() - PE0;
	double dotprod;
	for(int i=0; i<box.get_N(); i++){
		gsl_blas_ddot(vs[i],vs[i],&dotprod);
		dotprod=0.5*dotprod;
		nrg+=dotprod;
	}
	if(box.sym()>1){
		if(box.get_dim()==2){
			for(int i=0; i<box.get_N(); i++)
				nrg+=(0.5)*box.I(i)*wns[i]*wns[i];
		}
		else{//3D
			if(box.sym()==2){//I still number
				for(int i=0; i<box.get_N(); i++){
					gsl_blas_ddot(wvs[i],wvs[i],&dotprod);
					dotprod=0.5*box.I(i)*dotprod;
					nrg+=dotprod;
				}
			}
			else{//Imat
				//TODO do this for Imat
			}
		}
	}
	return nrg;
}

template<class P, class B>
inline void SingleMD<P, B>::next() {
	move_x();
	if(box.sym()>1){
		rot();
	}
	store_DU();
	if(box.sym()>1){
		if(box.get_dim()==2) store_Dtaun();
		else store_Dtauv();
	}
	(*pot).calc_U_F();
	move_v();
	if(box.sym()>1) move_w();
	//finally, store the new avg values and check if good
	double nrg = energy();
	avgU = (numstep*avgU+nrg)/((double)(numstep+1));
	avgU2 = (numstep*avgU2+nrg*nrg)/((double)(numstep+1));
	double stdev = std::sqrt(avgU2-avgU*avgU);
	if(stdev/avgU>std::pow(10,-4)){
		std::cout << "Energy is changing too much. Try a smaller step size.";
		std::exit(0);
	}
	//increase step number
	numstep++;
}

template<class P, class B>
inline void SingleMD<P, B>::move_x() {
	gsl_vector * delx = gsl_vector_alloc(box.get_dim());
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_memcpy(delx, (*pot).get_DU(i));
		gsl_vector_scale(delx,-0.5*delt);//(1/2)a*delt
		gsl_vector_add(delx,vs[i]);//v+(1/2)a*delt
		gsl_vector_scale(delx,delt);//v*delt+(1/2)a*delt^2; already in [0,1)
		gsl_vector_add(box.get_1_pos(i),delx);
		box.modify(i);
	}
	gsl_vector_free(delx);
}

template<class P, class B>
inline void SingleMD<P, B>::rot() {
	for(int i=0; i<box.get_N(); i++){
		if(box.get_dim()==2){
			double theta = wns[i]*delt-0.5*((*pot).get_Dtaun(i)/box.I(i))*delt*delt;
			double cose = std::cos(theta); double sine = std::sin(theta);
			gsl_matrix * rota = gsl_matrix_alloc(2,2);
			gsl_matrix_set(rota,0,0,cose);gsl_matrix_set(rota,0,1,-sine);
			gsl_matrix_set(rota,1,0,sine);gsl_matrix_set(rota,1,1,cose);
			gsl_vector * rotated = gsl_vector_alloc(2);
			gsl_blas_dgemv(CblasNoTrans,1.0,rota,box.get_u(i),0.0,rotated);
			gsl_vector_memcpy(box.get_u(i),rotated);
			gsl_vector_free(rotated);
			gsl_matrix_free(rota);
		}
		else{
			//TODO stuff
		}
		box.normalize(i);
	}
}

template<class P, class B>
inline void SingleMD<P, B>::store_DU() {
}

template<class P, class B>
inline void SingleMD<P, B>::store_Dtaun() {
}

template<class P, class B>
inline void SingleMD<P, B>::store_Dtauv() {
}

template<class P, class B>
inline void SingleMD<P, B>::move_v() {
}

template<class P, class B>
inline void SingleMD<P, B>::move_w() {
}

#endif /* SINGLEMD_H_ */
