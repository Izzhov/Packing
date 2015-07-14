/*
 * HarmPotNbrList.h
 *
 *  Created on: Jul 9, 2015
 *      Author: izzhov
 */

#ifndef HARMPOTNBRLIST_H_
#define HARMPOTNBRLIST_H_

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>

#include "Torus.h"
#include "mm.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <cmath>
#include <vector>

template <class B, class T>//B for Box, T for shape
class HarmPotNbrList {
private:
	double kk; double P;
	B& box;
	double U; std::vector<gsl_vector*> DU; //DU is for [0,1) coords
	std::vector<double> Dtaun; std::vector<gsl_vector*> Dtauv;
	double DL;
	std::vector<gsl_vector*> D_u; std::vector<gsl_vector*> D_v;

	std::vector<std::vector<int > > con;
	std::vector<T> skins;
	double L0;
public:
	HarmPotNbrList(B& bbox);

	void set_P(double P);//re-calculates energy etc after setting
	void set_DU(int i, int j, double x);
	void set_Dtaun(int i, double x);
	void set_Dtauv(int i, int j, double x);
	void set_D_u(int i, int j, double x);
	void set_D_v(int i, int j, double x);
	void set_DL(double newDL){DL = newDL;}

	double get_kk(){return kk;}
	double get_P(){return P;}
	double get_U(){return U;}
	gsl_vector * get_DU(int i){
		return DU[i];
	}
	gsl_vector * get_Dtauv(int i){
		return Dtauv[i];
	}
	gsl_vector * get_D_u(int i){
		return D_u[i];
	}
	gsl_vector * get_D_v(int i){
		return D_v[i];
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
	void calc_U_F(int i, int j, int k, int ncon);
	void calc_U_F(int i, int j, int ncon);
	double dUdl(double lijk, double Rij);

	bool skin_break(); //returns whether a particle is outside its skin
	void make_list(); //makes a new neighborlist+new skin particles
	int num_contacts(); //returns (twice) th'num o'contacts
	int num_skintacts(); //returns number that skins are said to be touchin'
	std::vector<int> get_cons(int i){return con[i];}//gives th'list of skintacts for ith

	double pressure();//calculates pressure
};

template<class B, class T>
inline HarmPotNbrList<B, T>::HarmPotNbrList(B& bbox):box(bbox){
	kk=1; P=std::pow(10,-4); //default values
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
	make_list();
	calc_U_F();
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::set_P(double P) {
	this->P = P;
	calc_U_F();
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::set_DU(int i, int j, double x) {
	gsl_vector_set(DU[i],j,x);
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::set_Dtaun(int i, double x) {
	Dtaun.at(i) = x;
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::set_Dtauv(int i, int j, double x) {
	gsl_vector_set(Dtauv[i],j,x);
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::set_D_u(int i, int j, double x) {
	gsl_vector_set(D_u[i],j,x);
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::set_D_v(int i, int j, double x) {
	gsl_vector_set(D_v[i],j,x);
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::calc_U_F() {
	if(skin_break()) make_list();
	bool usequadrants = (box.get_L()<=2*skins[box.get_lslmax()].lsl());
	double L = box.get_L(); int dim = box.get_dim();
	U = P*std::pow(L,dim); DL = dim*P*std::pow(L,dim-1);
	for(int i=0; i<box.get_N(); i++){
		set_Dtaun(i,0);
		for(int d=0; d<dim; d++){//initialize new forces to zero
			set_DU(i,d,0); set_Dtauv(i,d,0);
			set_D_u(i,d,0); set_D_v(i,d,0);
		}
		for(int j=0; j<con[i].size(); j++){
			if(usequadrants)
				for(int k=1; k<=mm::int_pow(2,dim);k++)
					for(int ncon=0; ncon<box.ncon(); ncon++)
						calc_U_F(i,con[i][j],k,ncon);
			else for(int ncon=0; ncon<box.ncon(); ncon++)
				calc_U_F(i,con[i][j],ncon);
		}
	}
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::calc_U_F(int i, int j, int k, int ncon) {
	if(i==j) return;
	double Rij = box.R(i,j);
	double lijk2 = box.ell2(i,j,k,ncon);
	if(lijk2>(Rij*Rij)) return;
	double lijk = sqrt(lijk2);
	if(i<j){//store energy
		U+=0.5*kk*(1-(lijk/Rij))*(1-(lijk/Rij));
	}
	double dudl = dUdl(lijk, Rij); double L = box.get_L();
	gsl_vector * f = box.ell_vec(i,j,k,ncon);//force direction (neg bc deriv)
	gsl_vector_scale(f,dudl/(L*lijk)); //force magnitude is good for [0,1)
	gsl_vector_add(DU[i],f);//adds this force to DU
	DL-=dudl*lijk/L;
	// TODO add torque rotations stuff
	if(box.sym()>1){
		gsl_vector * floc = box.F_loc(i,j,k,ncon);//location of force
		if(box.get_dim()==2){
			double newtaun=gsl_vector_get(floc,0)*gsl_vector_get(f,1)-
					gsl_vector_get(floc,1)*gsl_vector_get(f,0);//crossprod
			Dtaun.at(i)+=newtaun;
			gsl_vector * cloc = mm::cross(1,box.get_u(i));//for D_u (see 6-16-15 page 2)
			gsl_vector_scale(cloc,newtaun/box.I(i));
			gsl_vector_add(D_u[i],cloc);
			gsl_vector_free(cloc);
		}
		else if(box.sym()==2){//sym==2 but in 3d
			gsl_vector * newtauv = mm::cross(2,floc,f);
			gsl_vector_add(Dtauv[i],newtauv);
			gsl_vector * cloc = mm::cross(2,newtauv,box.get_u(i));
			gsl_vector_scale(cloc,1.0/box.I(i));
			gsl_vector_add(D_u[i],cloc);
			gsl_vector_free(newtauv); gsl_vector_free(cloc);
		}
		else{//sym==3 in 3d
			// TODO find new D_u and D_v with tauv and Imat
		}
		gsl_vector_free(floc);
	}
	gsl_vector_free(f);
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::calc_U_F(int i, int j, int ncon) {
	if(i==j) return;
	double Rij = box.R(i,j);
	double lijk2 = box.ell2(i,j,ncon);
	if(lijk2>(Rij*Rij)) return;
	double lijk = sqrt(lijk2);
	if(i<j){//store energy
		U+=0.5*kk*(1-(lijk/Rij))*(1-(lijk/Rij));
	}
	double dudl = dUdl(lijk, Rij); double L = box.get_L();
	gsl_vector * f = box.ell_vec(i,j,ncon);//force direction (neg bc deriv)
	gsl_vector_scale(f,dudl/(L*lijk)); //force magnitude is good for [0,1)
	gsl_vector_add(DU[i],f);//adds this force to DU
	DL-=dudl*lijk/L;
	// TODO add torque rotations stuff
	if(box.sym()>1){
		gsl_vector * floc = box.F_loc(i,j,ncon);//location of force
		if(box.get_dim()==2){
			double newtaun=gsl_vector_get(floc,0)*gsl_vector_get(f,1)-
					gsl_vector_get(floc,1)*gsl_vector_get(f,0);//crossprod
			Dtaun.at(i)+=newtaun;
			gsl_vector * cloc = mm::cross(1,box.get_u(i));//for D_u (see 6-16-15 page 2)
			gsl_vector_scale(cloc,newtaun/box.I(i));
			gsl_vector_add(D_u[i],cloc);
			gsl_vector_free(cloc);
		}
		else if(box.sym()==2){//sym==2 but in 3d
			gsl_vector * newtauv = mm::cross(2,floc,f);
			gsl_vector_add(Dtauv[i],newtauv);
			gsl_vector * cloc = mm::cross(2,newtauv,box.get_u(i));
			gsl_vector_scale(cloc,1.0/box.I(i));
			gsl_vector_add(D_u[i],cloc);
			gsl_vector_free(newtauv); gsl_vector_free(cloc);
		}
		else{//sym==3 in 3d
			// TODO find new D_u and D_v with tauv and Imat
		}
		gsl_vector_free(floc);
	}
	gsl_vector_free(f);
}

template<class B, class T>
inline double HarmPotNbrList<B, T>::dUdl(double lijk, double Rij) {
	return (kk/Rij)*(1-(lijk/Rij));
}

template<class B, class T>
inline bool HarmPotNbrList<B, T>::skin_break() {
	bool skinbreak=false;
	//first, update r's of part's
	double L=box.get_L();
	for(int i=0; i<box.get_N(); i++){
		skins[i].set_r(skins[i].get_r()*L/L0);
		if(box.get_shape(i).is_far(skins[i],L)){
			skinbreak = true; break;
		}
	}
	L0=L;//L0 is L from the previous time this was checked
	return skinbreak;
}

template<class B, class T>
inline void HarmPotNbrList<B, T>::make_list() {
	//std::cout << "making list" << std::endl;
	double L = box.get_L();
	//first need to free the gsl vectors of the skins
	for(int i=0; i<skins.size(); i++){
		gsl_vector_free(skins[i].get_pos());
		if(box.sym()>1){
			gsl_vector_free(skins[i].get_u());
			if(box.sym()==3) gsl_vector_free(skins[i].get_v());
		}
	}
	skins.clear(); con.clear(); L0=box.get_L();
	double ratio;//see 7-10-15 whiteback for ratio notes
	double avgcon=(double)num_contacts()/(double)box.get_N();
	double df;//degrees of freedom
	if(box.get_dim()==2) df=box.get_dim()+box.sym()-1;
	else df=box.get_dim()+box.sym()-(box.sym()==1);
	if(avgcon<1) ratio=2.0;
	else if(avgcon>df+1) ratio=1.05;
	else ratio=(1.05*(avgcon-1)+2*(df+1-avgcon))/df;
	if(ratio<1){
		std::cout << "ratio < 1" << std::endl;
		exit(0);
	}
	//make skins
	for(int i=0; i<box.get_N(); i++){
		T nextshape(box.get_shape_ptr(i),ratio);
		skins.push_back(nextshape);
	}
	//make con
	bool usequadrants = (box.get_L()<=2*skins[box.get_lslmax()].lsl());
	for(int i=0; i<box.get_N(); i++){
		std::vector<int> nex; //list of contacts for i
		con.push_back(nex); //nex will be at i
		for(int j=0; j<box.get_N(); j++){
			if(i==j){ //can't touch itself
				if(i==(box.get_N()-1)) break;//this means we're done
				else j++;//move on to the next one if not done
			}
			if(usequadrants){
				for(int k=1; k<=mm::int_pow(2,box.get_dim());k++)
					for(int ncon=0; ncon<box.ncon(); ncon++)
						if(L*std::sqrt(skins.at(i).ell2(skins.at(j),k,L,ncon))<=
								skins.at(i).get_r()+skins.at(j).get_r()){
							if(con[i].size()>0){
								if(con.at(i).at(con.at(i).size()-1)!=j){
									//std::cout << con.at(i).at(con.at(i).size()-1) << std::endl;
									con.at(i).push_back(j);
								}
							}
							else{
								con.at(i).push_back(j);
							}
						}
			}
			else for(int ncon=0; ncon<box.ncon(); ncon++){
				if(L*std::sqrt(skins.at(i).ell2(skins.at(j),L,ncon))<=
						skins.at(i).get_r()+skins.at(j).get_r()){
					if(con[i].size()>0){
						if(con.at(i).at(con.at(i).size()-1)!=j){
							//std::cout << con.at(i).at(con.at(i).size()-1) << std::endl;
							con.at(i).push_back(j);
						}
					}
					else{
						con.at(i).push_back(j);
					}
				}
			}
		}
	}
	/*
	//check if number of con-tacts is >= num_contacts (it should be)
	double con_tacts = 0;
	for(int i=0; i<box.get_N(); i++) con_tacts+=con[i].size();
	if(con_tacts<num_contacts()){
		std::cout << "invisible neighbors" << std::endl;
		exit(0);
	}
	*/
}

template<class B, class T>
inline int HarmPotNbrList<B, T>::num_contacts() {
	bool usequadrants = (box.get_L()<=2*box.lsl(box.get_lslmax()));
	if(usequadrants) std::cout << "we're using quadrants now?" << std::endl;
	int numcontacts=0;
	for(int i=0; i<box.get_N(); i++){
		for(int j=0; j<box.get_N(); j++){
			if(i==j){ //can't touch itself
				if(i==(box.get_N()-1)) break;//this means we're done
				else j++;//move on to the next one if not done
			}
			if(usequadrants){
				for(int k=1; k<=mm::int_pow(2,box.get_dim());k++)
					for(int ncon=0; ncon<box.ncon(); ncon++)
						if(std::sqrt(box.ell2(i,j,k,ncon))<=box.R(i,j))
							numcontacts++;
			}
			else for(int ncon=0; ncon<box.ncon(); ncon++){
				if(std::sqrt(box.ell2(i,j,ncon))<=box.R(i,j)){
					numcontacts++;
				}
			}
		}
	}
	return numcontacts;
}

template<class B, class T>
inline int HarmPotNbrList<B, T>::num_skintacts(){
	int skintacts = 0;
	for(int i=0; i<box.get_N(); i++) skintacts+=con[i].size();
	return skintacts;
}

template<class B, class T>
inline double HarmPotNbrList<B, T>::pressure() {
	double L = box.get_L(); int dim = box.get_dim();
	return P-DL/(dim*std::pow(L,dim-1));
}

#endif /* HARMPOTNBRLIST_H_ */
