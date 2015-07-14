/*
 * Sphere.cpp
 *
 *  Created on: Jun 8, 2015
 *      Author: izzhov
 */

#include "Sphere.h"

Sphere::Sphere(){}

Sphere::Sphere(gsl_vector* pos, std::vector<gsl_vector*> size) {
	set_pos(pos); set_r(gsl_vector_get(size[0],0));
}

Sphere::Sphere(gsl_vector * pos, gsl_vector * u,
		std::vector<gsl_vector*> size){}

Sphere::Sphere(gsl_vector * pos, gsl_vector * u,
		gsl_vector * v, std::vector<gsl_vector*> size){}

void Sphere::set_pos(gsl_vector* pos) {
	this->pos = pos;
}

void Sphere::set_r(double r) {
	this->r = r;
}

void Sphere::set_pos_coord(int i, double x){
	gsl_vector_set(pos,i,mm::pmod(x,1));
}

gsl_vector* Sphere::get_pos() {
	return pos;
}

double Sphere::get_r() {
	return r;
}

double Sphere::get_pos_coord(int i){
	return gsl_vector_get(pos,i);
}

double Sphere::volume() {
	if(pos->size == 2) return M_PI*r*r;
	else return (4.0/3.0)*M_PI*r*r*r;
}

double Sphere::lsl(){
	return 2*r;
}

gsl_vector * Sphere::ell_vec(Sphere s, int k, double L, int ncon){
	return mm::rel(pos,s.get_pos(),k);
}

double Sphere::ell2(Sphere s, int k, double L, int ncon) {
	double dd;//the distance
	gsl_vector * reld = ell_vec(s,k,L,ncon);
	gsl_blas_ddot(reld,reld,&dd);
	gsl_vector_free(reld);
	return dd;
}

gsl_vector * Sphere::ell_vec(Sphere s, double L, int ncon){
	return mm::rel(pos,s.get_pos());
}

double Sphere::ell2(Sphere s, double L, int ncon) {
	double dd;//the distance
	gsl_vector * reld = ell_vec(s,L,ncon);
	gsl_blas_ddot(reld,reld,&dd);
	gsl_vector_free(reld);
	return dd;
}

bool Sphere::touch(Sphere s, double L){
	bool dotheytouch = false;
	for(int k=1; k<=mm::int_pow(2,pos->size);k++){
		if(std::sqrt(ell2(s,k,L,0))<=s.get_r()+get_r()) dotheytouch = true;
	}
	return dotheytouch;
}

void Sphere::normalize() {
}

double Sphere::max_d() {
}

gsl_vector* Sphere::get_u() {
}

gsl_vector* Sphere::get_v() {
}

double Sphere::I() {
}

gsl_vector * Sphere::F_loc(Sphere s, int k, double L, int ncon) {
	gsl_vector * floc = gsl_vector_calloc(pos->size);
	return floc;
}

gsl_vector * Sphere::F_loc(Sphere s, double L, int ncon) {
	gsl_vector * floc = gsl_vector_calloc(pos->size);
	return floc;
}
