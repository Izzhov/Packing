/*
 * MultiMD.h
 *
 *  Created on: Jun 25, 2015
 *      Author: izzhov
 */

#ifndef MULTIMD_H_
#define MULTIMD_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include "mm.h"
#include "HarmPot.h"
#include "Torus.h"
#include "TimeRNGNormal.h"
#include "SingleMD.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <vector>
#include <cmath>

template <class P, class B>//P is potential, B is box
class MultiMD {
private:
	SingleMD<P,B>& smd;
	std::vector<SingleMD<P,B>*> MDs;
	B box;
	double autocorr;
	std::vector<std::vector<gsl_vector*> > initvs;
	std::vector<std::vector<double> > initwns;
	std::vector<std::vector<gsl_vector*> > initwvs;
	std::ofstream dstream;
	double dDt;//1/lgfreq
	double Dt;//2pi/smfreq
	int numdels;//number of dDt's in Dt (aka INT)
	//# of elements in MDs aka # of dotprods to avg over in autocorr
	int autoavgnum;
	int dDt_o_delt;//number of delts in dDt
public:
	MultiMD(SingleMD<P,B>& ssmd, double smfreq, double lgfreq);

	void populate();

	void store_vs();
	void store_wns();
	void store_wvs();

	void find_autocorr();
	void add_autocorr(int i);
};

template<class P, class B>
inline MultiMD<P, B>::MultiMD(SingleMD<P, B>& ssmd, double smfreq,
		double lgfreq):smd(ssmd) {
	autoavgnum = 1000; dDt_o_delt = 100;
	box = smd.get_box(); autocorr = 0;
	dDt = 1.0/lgfreq; Dt = 2.0*M_PI/smfreq;
	numdels = std::floor(Dt/dDt)+1;
	smd.set_delt(dDt/(double)dDt_o_delt);
	populate();
}

template<class P, class B>
inline void MultiMD<P, B>::populate() {
	for(int i=0; i<autoavgnum; i++){
		//store the veloc's
		store_vs();
		if(box.sym()>1){
			if(box.get_dim()==2) store_wns();
			else store_wvs();
		}
		//store the MD state
		SingleMD<P,B> * nextsmd = new SingleMD<P,B>(smd,0);
		MDs.push_back(nextsmd);
		//iterate to next one (unless final storage)
		if(i<autoavgnum-1) for(int j=0; j<dDt_o_delt*numdels; j++)
			smd.next();
	}
}

template<class P, class B>
inline void MultiMD<P, B>::store_vs() {
	std::vector<gsl_vector*> nextvs;
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextv = gsl_vector_alloc(box.get_dim());
		gsl_vector_memcpy(nextv,smd.get_vs(i));
		nextvs.push_back(nextv);
	}
	initvs.push_back(nextvs);
}

template<class P, class B>
inline void MultiMD<P, B>::store_wns() {
	std::vector<double> nextwns;
	for(int i=0; i<box.get_N(); i++) nextwns.push_back(smd.get_wns(i));
	initwns.push_back(nextwns);
}

template<class P, class B>
inline void MultiMD<P, B>::store_wvs() {
	std::vector<gsl_vector*> nextwvs;
	for(int i=0; i<box.get_N(); i++){
		gsl_vector * nextwv = gsl_vector_alloc(box.get_dim());
		gsl_vector_memcpy(nextwv,smd.get_wvs(i));
		nextwvs.push_back(nextwv);
	}
	initwvs.push_back(nextwvs);
}

template<class P, class B>
inline void MultiMD<P, B>::find_autocorr() {
	remove("autocorr.txt");
	dstream.open("autocorr.txt");
	dstream << std::setprecision(16);
	for(int i=1; i<dDt_o_delt*numdels+1; i++){
		dstream << i*smd.get_delt() << " ";
		for(int j=0; j<autoavgnum; j++){
			(*MDs.at(j)).next();
			add_autocorr(j);
		}
		autocorr = autocorr/(double)(box.get_N()*autoavgnum);
		dstream << autocorr << std::endl;
		autocorr = 0;
	}
	dstream.close();
}

template<class P, class B>
inline void MultiMD<P, B>::add_autocorr(int i) {
	double dotprod;
	for(int j=0; j<box.get_N(); j++){
		gsl_blas_ddot((*MDs.at(i)).get_vs(j),initvs.at(i).at(j),&dotprod);
		autocorr+=dotprod;
		if(box.sym()>1){
			if(box.get_dim()==2){
				autocorr+=box.I(j)*(*MDs.at(i)).get_wns(j)*initwns.at(i).at(j);
			}
			else{
				if(box.sym()==2){
					gsl_blas_ddot((*MDs.at(i)).get_wvs(j),initwvs.at(i).at(j),&dotprod);
					dotprod*=box.I(j);
					autocorr+=dotprod;
				}
				else{
					//TODO use Imat
				}
			}
		}
	}
}

#endif /* MULTIMD_H_ */
