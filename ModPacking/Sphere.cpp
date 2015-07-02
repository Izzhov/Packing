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

gsl_vector * Sphere::ell_vec(Sphere s, int k, double L){
	return mm::rel(pos,s.get_pos(),k);
}

double Sphere::ell2(Sphere s, int k, double L) {
	double dd;//the distance
	gsl_vector * reld = ell_vec(s,k,L);
	gsl_blas_ddot(reld,reld,&dd);
	gsl_vector_free(reld);
	return dd;
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

gsl_vector * Sphere::F_loc(Sphere s, int k, double L) {
	gsl_vector * floc = gsl_vector_calloc(pos->size);
	return floc;
}

