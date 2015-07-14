/*
 * Spherocyl.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: izzhov
 */

#include "Spherocyl.h"

Spherocyl::Spherocyl() {}

Spherocyl::Spherocyl(gsl_vector* pos, std::vector<gsl_vector*> size) {}

Spherocyl::Spherocyl(gsl_vector* pos, gsl_vector* u,
		std::vector<gsl_vector*> size) {
	set_pos(pos); set_u(u); set_r(gsl_vector_get(size[0],0));
	set_a(gsl_vector_get(size[1],0));
}

Spherocyl::Spherocyl(gsl_vector* pos, gsl_vector* u, gsl_vector* v,
		std::vector<gsl_vector*> size) {}

void Spherocyl::set_pos(gsl_vector* pos) {
	this->pos = pos;
}

void Spherocyl::set_u(gsl_vector* u) {
	this->u = u;
}

void Spherocyl::set_r(double r) {
	this->r = r;
}

void Spherocyl::set_a(double a) {
	this->a = a;
}

void Spherocyl::set_pos_coord(int i, double x) {
	gsl_vector_set(pos,i,mm::pmod(x,1));
}

void Spherocyl::set_u_coord(int i, double x) {
	gsl_vector_set(u,i,x);
}

double Spherocyl::sphere_volume(){
	if(pos->size == 2) return (M_PI*r*r);
	else return (4.0/3.0)*M_PI*r*r*r;
}

double Spherocyl::cyl_volume(){
	if(pos->size == 2) return 4*(a-1)*r*r;
	else return M_PI*r*r*(2*(a-1)*r);
}

double Spherocyl::volume() {
	return sphere_volume()+cyl_volume();
}

double Spherocyl::I() {
	double V = volume(); double S = sphere_volume(); double C = cyl_volume();
	double mS = S/V; double mC = C/V;
	if(pos->size == 2){
		double h = 2*(a-1)*r; double w = 2*r;
		double IS = (r*r/2.0+h*h/4.0);
		double IC = (h*h+w*w)/12.0;
		return mS*IS+mC*IC;
	}
	else{//3D
		double L = 2*(a-1)*r;
		double IS = ((2.0/5.0)*r*r+L*L/4.0);
		double IC = (1.0/4.0)*r*r+(1.0/12.0)*L*L;
		return mS*IS+mC*IC;
	}
}

double Spherocyl::max_d() {
	return (a-1)*r;
}

double Spherocyl::lsl(){
	return 2*a*r;
}

gsl_vector * Spherocyl::lisljs(Spherocyl s, int k, double L){
	gsl_vector * thels = gsl_vector_calloc(2);
	double lis = 0; double ljs = 0;
	double li = 2*max_d(); double lj = 2*s.max_d();
	double uiuj; gsl_blas_ddot(get_u(),s.get_u(),&uiuj);
	gsl_vector * rij = mm::rel(pos,s.get_pos(),k);
	gsl_vector_scale(rij,L);
	double uirij; gsl_blas_ddot(get_u(), rij, &uirij);
	double ujrij; gsl_blas_ddot(s.get_u(),rij,&ujrij);
	if((1-uiuj*uiuj)<std::pow(10,-8)){//then spherocyls are parallel
		if(std::abs(uirij)>(li+lj)/2.0){//force only applied in one place
			lis = mm::sgn(uirij)*li/2.0; ljs = mm::sgn(-ujrij)*lj/2.0;
		}
		else{//force applied in middle of 2 forcends (see6-16-15 page 2 for algo)
			//see 7-1-15 for correction (ujrij -> -ujrij)
			double posi = std::min(li/2.0,uirij+lj/2.0);
			double negi = std::max(-li/2.0,uirij-lj/2.0);
			double posj = std::min(lj/2.0,-ujrij+li/2.0);
			double negj = std::max(-lj/2.0,-ujrij-li/2.0);
			lis = (posi+negi)/2.0; ljs = (posj+negj)/2.0;
		}
	}
	else{//not parallel, use regular algorithm from paper
		double lip = (uirij-uiuj*ujrij)/(1-uiuj*uiuj);
		double ljp = (uiuj*uirij-ujrij)/(1-uiuj*uiuj);
		double bli = std::abs(lip)-li/2.0;
		double blj = std::abs(ljp)-lj/2.0;
		if(bli<=0 && blj <=0){
			lis = lip; ljs = ljp;
		}
		else{
			if(blj>=bli){
				ljs = mm::sign(lj/2.0,ljp);
				lis = std::max(-li/2.0,std::min(uirij+ljs*uiuj,li/2.0));
			}
			else{
				lis = mm::sign(li/2.0,lip);
				ljs = std::max(-lj/2.0,std::min(-ujrij+lis*uiuj,lj/2.0));
			}
		}
	}
	gsl_vector_set(thels,0,lis); gsl_vector_set(thels,1,ljs);
	gsl_vector_free(rij);
	return thels;
}

gsl_vector* Spherocyl::F_loc(Spherocyl s, int k, double L, int ncon) {
	gsl_vector * floc = gsl_vector_alloc(pos->size);
	gsl_vector * thels = lisljs(s,k,L);
	double li = gsl_vector_get(thels,0);
	gsl_vector_memcpy(floc,get_u()); gsl_vector_scale(floc,li);
	gsl_vector_free(thels);
	return floc;
}

gsl_vector* Spherocyl::ell_vec(Spherocyl s, int k, double L, int ncon) {
	gsl_vector * rij = mm::rel(pos,s.get_pos(),k);
	gsl_vector * uu = gsl_vector_alloc(rij->size);
	gsl_vector * thels = lisljs(s,k,L);
	double li = gsl_vector_get(thels,0); double lj = gsl_vector_get(thels,1);
	gsl_vector_memcpy(uu,get_u()); gsl_vector_scale(uu,-li/L);
	gsl_vector_add(rij,uu);
	gsl_vector_memcpy(uu,s.get_u()); gsl_vector_scale(uu,lj/L);
	gsl_vector_add(rij,uu);
	gsl_vector_free(uu); gsl_vector_free(thels);
	return rij;
}

double Spherocyl::ell2(Spherocyl s, int k, double L, int ncon) {
	double dd;//the distance
	gsl_vector * reld = ell_vec(s,k,L, ncon);
	gsl_blas_ddot(reld,reld,&dd);
	gsl_vector_free(reld);
	return dd;
}

gsl_vector * Spherocyl::lisljs(Spherocyl s, double L){
	gsl_vector * thels = gsl_vector_calloc(2);
	double lis = 0; double ljs = 0;
	double li = 2*max_d(); double lj = 2*s.max_d();
	double uiuj; gsl_blas_ddot(get_u(),s.get_u(),&uiuj);
	gsl_vector * rij = mm::rel(pos,s.get_pos());
	gsl_vector_scale(rij,L);
	double uirij; gsl_blas_ddot(get_u(), rij, &uirij);
	double ujrij; gsl_blas_ddot(s.get_u(),rij,&ujrij);
	if((1-uiuj*uiuj)<std::pow(10,-8)){//then spherocyls are parallel
		if(std::abs(uirij)>(li+lj)/2.0){//force only applied in one place
			lis = mm::sgn(uirij)*li/2.0; ljs = mm::sgn(-ujrij)*lj/2.0;
		}
		else{//force applied in middle of 2 forcends (see6-16-15 page 2 for algo)
			//see 7-1-15 for correction (ujrij -> -ujrij)
			double posi = std::min(li/2.0,uirij+lj/2.0);
			double negi = std::max(-li/2.0,uirij-lj/2.0);
			double posj = std::min(lj/2.0,-ujrij+li/2.0);
			double negj = std::max(-lj/2.0,-ujrij-li/2.0);
			lis = (posi+negi)/2.0; ljs = (posj+negj)/2.0;
		}
	}
	else{//not parallel, use regular algorithm from paper
		double lip = (uirij-uiuj*ujrij)/(1-uiuj*uiuj);
		double ljp = (uiuj*uirij-ujrij)/(1-uiuj*uiuj);
		double bli = std::abs(lip)-li/2.0;
		double blj = std::abs(ljp)-lj/2.0;
		if(bli<=0 && blj <=0){
			lis = lip; ljs = ljp;
		}
		else{
			if(blj>=bli){
				ljs = mm::sign(lj/2.0,ljp);
				lis = std::max(-li/2.0,std::min(uirij+ljs*uiuj,li/2.0));
			}
			else{
				lis = mm::sign(li/2.0,lip);
				ljs = std::max(-lj/2.0,std::min(-ujrij+lis*uiuj,lj/2.0));
			}
		}
	}
	gsl_vector_set(thels,0,lis); gsl_vector_set(thels,1,ljs);
	gsl_vector_free(rij);
	return thels;
}

gsl_vector* Spherocyl::F_loc(Spherocyl s, double L, int ncon) {
	gsl_vector * floc = gsl_vector_alloc(pos->size);
	gsl_vector * thels = lisljs(s,L);
	double li = gsl_vector_get(thels,0);
	gsl_vector_memcpy(floc,get_u()); gsl_vector_scale(floc,li);
	gsl_vector_free(thels);
	return floc;
}

gsl_vector* Spherocyl::ell_vec(Spherocyl s, double L, int ncon) {
	gsl_vector * rij = mm::rel(pos,s.get_pos());
	gsl_vector * uu = gsl_vector_alloc(rij->size);
	gsl_vector * thels = lisljs(s,L);
	double li = gsl_vector_get(thels,0); double lj = gsl_vector_get(thels,1);
	gsl_vector_memcpy(uu,get_u()); gsl_vector_scale(uu,-li/L);
	gsl_vector_add(rij,uu);
	gsl_vector_memcpy(uu,s.get_u()); gsl_vector_scale(uu,lj/L);
	gsl_vector_add(rij,uu);
	gsl_vector_free(uu); gsl_vector_free(thels);
	return rij;
}

double Spherocyl::ell2(Spherocyl s, double L, int ncon) {
	double dd;//the distance
	gsl_vector * reld = ell_vec(s,L, ncon);
	gsl_blas_ddot(reld,reld,&dd);
	gsl_vector_free(reld);
	return dd;
}

bool Spherocyl::touch(Spherocyl s, double L){
	bool dotheytouch = false;
	for(int k=1; k<=mm::int_pow(2,pos->size);k++){
		if(std::sqrt(ell2(s,k,L,0))<=s.get_r()+get_r()) dotheytouch = true;
	}
	return dotheytouch;
}

void Spherocyl::normalize() {
	double unorm = gsl_blas_dnrm2(u);
	gsl_vector_scale(u,1.0/unorm);
}

gsl_vector* Spherocyl::get_v() {
}
