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
	TimeRNGNormal trng;
	double PE0; //minimal PE
	double avgU; double avgU2;//for finding avg and std-dev
	double skT; //10^(-8)
	B box; P * pot;
	std::vector<gsl_vector*> oldDU;
	std::vector<double> oldDtaun; std::vector<gsl_vector*> oldDtauv;
	std::vector<gsl_vector*> vs;
	std::vector<double> wns; std::vector<gsl_vector*> wvs;
	double delt; int numstep;
public:
	SingleMD(B bbox, double ddelt, TimeRNGNormal ttrng);
	//copy constructor, ndef is so not mixed up with default copy
	SingleMD(SingleMD<P,B>& copy, int ndef);
	//set velocities (don't use when sym>1)
	SingleMD(B bbox, double ddelt,
			std::vector<gsl_vector*> vvs, TimeRNGNormal ttrng);
	//set velocities based on an eigenvector (don't use when sym>1)
	SingleMD(B bbox, double ddelt, std::vector<double> eig);

	void set_delt(double delt){this->delt = delt;}

	double get_PE0(){return PE0;}
	double get_avgU(){return avgU;}
	double get_avgU2(){return avgU2;}
	double get_skT(){return skT;}
	B get_box(){return box;}
	B& get_boxaddress(){return box;}
	gsl_vector* get_oldDU(int i){return oldDU[i];}
	std::vector<double> get_oldDtaun(){return oldDtaun;}
	gsl_vector* get_oldDtauv(int i){return oldDtauv[i];}
	gsl_vector* get_vs(int i){return vs[i];}
	std::vector<double> get_wns(){return wns;}
	gsl_vector* get_wvs(int i){return wvs[i];}
	double get_delt(){return delt;}
	int get_numstep(){return numstep;}
	P get_pot(){return (*pot);}
	TimeRNGNormal get_trng(){return trng;}

	double get_wns(int i){return wns[i];}

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
inline SingleMD<P, B>::SingleMD(B bbox, double ddelt,
		TimeRNGNormal ttrng){
	trng = ttrng;
	box = bbox; delt = ddelt; numstep=1;
	pot = new P(box);
	PE0 = (*pot).get_U(); skT = std::pow(10,-6);
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
inline SingleMD<P, B>::SingleMD(SingleMD<P,B>& copy,
		int ndef){
	PE0 = copy.get_PE0(); avgU = copy.get_avgU(); avgU2 = copy.get_avgU2();
	skT = copy.get_skT(); oldDtaun = copy.get_oldDtaun(); wns = copy.get_wns();
	delt = copy.get_delt(); numstep = copy.get_numstep();
	B bbox(copy.get_boxaddress(),0); box = bbox;
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
inline SingleMD<P, B>::SingleMD(B bbox, double ddelt,
		std::vector<gsl_vector*> vvs, TimeRNGNormal ttrng){
	trng=ttrng;
	box = bbox; delt = ddelt; numstep=1;
	pot = new P(box);
	PE0 = (*pot).get_U(); skT = std::pow(10,-6);
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextDU = gsl_vector_alloc(box.get_dim());
		gsl_vector * nextDtauv = gsl_vector_alloc(box.get_dim());
		gsl_vector_set_zero(nextDU); gsl_vector_set_zero(nextDtauv);
		oldDU.push_back(nextDU); oldDtauv.push_back(nextDtauv);
		oldDtaun.push_back(0);
	}
	vs = vvs;
	double nrg = energy();
	avgU = nrg; avgU2 = nrg*nrg;
}

template<class P, class B>
inline SingleMD<P, B>::SingleMD(B bbox, double ddelt, std::vector<double> eig){
	box = bbox; delt = ddelt; numstep=1;
	pot = new P(box);
	PE0 = (*pot).get_U(); skT = std::pow(10,-6);
	//store the forces
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextDU = gsl_vector_alloc(box.get_dim());
		gsl_vector * nextDtauv = gsl_vector_alloc(box.get_dim());
		gsl_vector_set_zero(nextDU); gsl_vector_set_zero(nextDtauv);
		oldDU.push_back(nextDU); oldDtauv.push_back(nextDtauv);
		oldDtaun.push_back(0);
	}
	//store the velocities
	int dim = box.get_dim();
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextv = gsl_vector_alloc(dim);
		for(int d=0; d<dim; d++){
			gsl_vector_set(nextv,d,skT*eig[i*dim+d]);
		}
		vs.push_back(nextv);
	}
	double nrg = energy();
	avgU = nrg; avgU2 = nrg*nrg;
}

template<class P, class B>
inline void SingleMD<P, B>::populate_v(){
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextv = gsl_vector_alloc(box.get_dim());
		for(int d=0; d<box.get_dim(); d++)
			gsl_vector_set(nextv,d,skT*trng.num());
		vs.push_back(nextv);
	}
}

template<class P, class B>
inline void SingleMD<P, B>::populate_wn(){
	for(int i=0; i<box.get_N(); i++){
		double sI = std::sqrt(box.I(i));
		wns.push_back((skT/sI)*trng.num());
	}
}

template<class P, class B>
inline void SingleMD<P, B>::populate_wv() {
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
		//need to convert to real energy units
		dotprod=0.5*dotprod*box.get_L()*box.get_L();
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
	if(stdev/(avgU+PE0)>std::pow(10,-4)){
		std::cout << "Energy is changing too much. Try a smaller step size." << std::endl;
		std::cout << "AvgU+PE0: " << avgU+PE0 << std::endl;
		std::cout << "StdDev: " << stdev << std::endl;
		std::cout << "NewNRG: " << energy() << std::endl;
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
	gsl_vector * delth = gsl_vector_alloc(box.get_dim());
	gsl_vector * cross = gsl_vector_alloc(box.get_dim());
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
		else{//dim==3
			if(box.sym()==2){//only u
				gsl_vector_memcpy(delth,(*pot).get_Dtauv(i));
				gsl_vector_scale(delth,-0.5*delt/box.I(i));//(1/2)a*delt
				gsl_vector_add(delth,wvs[i]);//v+(1/2)a*delt
				gsl_vector_scale(delth,delt);//v*delt+(1/2)a*delt^2
				double theta = gsl_blas_dnrm2(delth);
				double sine = std::sin(theta); double cose = std::cos(theta);
				gsl_vector_scale(delth,(1.0)/gsl_blas_dnrm2(delth));
				cross = mm::cross(2,delth,box.get_u(i));
				double dotprod; gsl_blas_ddot(delth,box.get_u(i),&dotprod);
				gsl_vector_scale(delth,(1-cose)*dotprod);
				gsl_vector_scale(cross,sine);
				gsl_vector_scale(box.get_u(i),cose);
				gsl_vector_add(box.get_u(i),cross);
				gsl_vector_add(box.get_u(i),delth);
			}
			else{
				//TODO when box.sym()==3
			}
		}
		box.normalize(i);
	}
	gsl_vector_free(delth);
	gsl_vector_free(cross);
}

template<class P, class B>
inline void SingleMD<P, B>::store_DU() {
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_memcpy(oldDU[i],(*pot).get_DU(i));
	}
}

template<class P, class B>
inline void SingleMD<P, B>::store_Dtaun() {
	for(int i=0; i<box.get_N(); i++){
		oldDtaun[i] = (*pot).get_Dtaun(i);
	}
}

template<class P, class B>
inline void SingleMD<P, B>::store_Dtauv() {
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_memcpy(oldDtauv[i],(*pot).get_Dtauv(i));
	}
}

template<class P, class B>
inline void SingleMD<P, B>::move_v() {
	gsl_vector * aold = gsl_vector_alloc(box.get_dim());
	gsl_vector * anew = gsl_vector_alloc(box.get_dim());
	for(int i=0; i<box.get_N(); i++){
		gsl_vector_memcpy(anew, (*pot).get_DU(i));
		gsl_vector_scale(anew,-0.5*delt);//(1/2)anew*delt
		gsl_vector_memcpy(aold,oldDU[i]);
		gsl_vector_scale(aold,-0.5*delt);//(1/2)aold*delt
		gsl_vector_add(vs[i],anew);
		gsl_vector_add(vs[i],aold);
	}
	gsl_vector_free(aold);
	gsl_vector_free(anew);
}

template<class P, class B>
inline void SingleMD<P, B>::move_w() {
	gsl_vector * aold = gsl_vector_alloc(box.get_dim());
	gsl_vector * anew = gsl_vector_alloc(box.get_dim());
	for(int i=0; i<box.get_N(); i++){
		if(box.get_dim()==2){
			wns[i]-=0.5*(oldDtaun[i]+(*pot).get_Dtaun(i))*(delt/box.I(i));
		}
		else{
			if(box.sym()==2){
				gsl_vector_memcpy(anew,(*pot).get_Dtauv(i));
				gsl_vector_scale(anew,-0.5*delt/box.I(i));
				gsl_vector_memcpy(aold,oldDtauv[i]);
				gsl_vector_scale(aold,-0.5*delt/box.I(i));
				gsl_vector_add(wvs[i],anew);
				gsl_vector_add(wvs[i],aold);
			}
			else{
				//TODO when Imat
			}
		}
	}
	gsl_vector_free(aold);
	gsl_vector_free(anew);
}

#endif /* SINGLEMD_H_ */
