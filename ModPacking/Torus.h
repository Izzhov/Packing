/*
 * Torus.h
 *
 *  Created on: Jun 8, 2015
 *      Author: izzhov
 */

#ifndef TORUS_H_
#define TORUS_H_

#include "Sphere.h"
#include "Sizes.h"
#include "TimeRNG01.h"
#include "mm.h"

#include <gsl/gsl_vector.h>

#include <vector>
#include <cmath>

template <class T>
class Torus {
private:
	int N;
	double L;
	Sizes ssizes;
	TimeRNG01 trng;
	std::vector<T> shapes;
	double partvol;
	int dim;//dimension of torus
public:
	Torus();
	Torus(int N, std::vector<std::vector<gsl_vector*> > szs,
			std::vector<double> wts, int dim, double phi);
	//to set seed manually
	Torus(int N, std::vector<std::vector<gsl_vector*> > szs,
			std::vector<double> wts, int dim, double phi, int seed);
	//to populate with pre-set locations and sizes
	Torus(int N, double L, std::vector<std::vector<gsl_vector*> > locsusvs,
			std::vector<std::vector<gsl_vector*> > fullsizes, int dim);
	//copy constructor, ndef is because I don't wanna confuse this w the default
	Torus(Torus<T>& copy, int ndef);

	int sym(){return T::sym(dim);}

	int get_seed(){return trng.get_seed();}

	//populates randomly
	void populate();

	void find_set_L(double phi);

	void set_N(int N);
	void set_L(double L);
	void set_dim(int dim);
	//particle # i, coordinate index j, is set to x (note: sets the [0,1) coord)
	void set_pos_coord(int i, int j, double x);
	void set_u_coord(int i, int j, double x);
	void set_v_coord(int i, int j, double x);

	int get_N(){return N;}
	double get_L(){return L;}
	int get_dim(){return dim;}
	double get_partvol(){return partvol;}
	Sizes get_sizes(){return ssizes;}
	TimeRNG01 get_trng(){return trng;}

	double get_r(int i){return shapes[i].get_r();}
	double get_max_d(int i){return shapes[i].max_d();}

	gsl_vector * get_1_pos(int i){return shapes[i].get_pos();}
	gsl_vector * get_u(int i){return shapes[i].get_u();}
	gsl_vector * get_v(int i){return shapes[i].get_v();}
	double get_1_pos_coord(int i, int j){return shapes[i].get_pos_coord(j);}
	double get_pos_coord(int i, int j){return L*get_1_pos_coord(i,j);}
	double get_u_coord(int i, int j){return shapes[i].get_u_coord(j);}
	double get_v_coord(int i, int j){return shapes[i].get_v_coord(j);}

	T get_shape(int i){return shapes[i];}

	double R(int i, int j);//distance within which they're overlapping
	double I(int i);//moment of inertia

	//location r-vec of the force vector
	gsl_vector * F_loc(int i, int j, int k);
	//returns in [0,1) space
	gsl_vector * ell_1_vec(int i, int j, int k);
	//returns in [0,L) space, is the direction of the force vector
	gsl_vector * ell_vec(int i, int j, int k);
	double ell2(int i, int j, int k);

	//packing faction
	double pack_frac();

	//normalizes the i-th shape
	void normalize(int i);

	//puts the coord back in [0,1)
	void modify(int i);
};

template<class T>
Torus<T>::Torus() {
	N=0;L=1;partvol=0;dim=1;
}

template<class T>
Torus<T>::Torus(int N, std::vector<std::vector<gsl_vector*> > szs,
		std::vector<double> wts, int dim, double phi) {
	partvol=0; set_N(N); set_dim(dim);
	Sizes newssizes(N, szs, wts); ssizes = newssizes;
	populate();
	for(int i=0; i<N; i++) partvol += shapes[i].volume();
	find_set_L(phi);
}

template<class T>
Torus<T>::Torus(int N, std::vector<std::vector<gsl_vector*> > szs,
		std::vector<double> wts, int dim, double phi, int seed) {
	TimeRNG01 newtrng(seed); trng = newtrng;
	partvol=0; set_N(N); set_dim(dim);
	Sizes newssizes(N, szs, wts); ssizes = newssizes;
	populate();
	for(int i=0; i<N; i++) partvol += shapes[i].volume();
	find_set_L(phi);
}

template<class T>
Torus<T>::Torus(int N, double L, std::vector<std::vector<gsl_vector*> > locsusvs,
		std::vector<std::vector<gsl_vector*> > fullsizes, int dim){
	partvol=0; set_N(N); set_dim(dim); set_L(L);
	for(int i=0; i<N; i++){
		T nextshape;
		if(sym()==1){
			T nextshape2(locsusvs[i][0], fullsizes[i]);
			nextshape = nextshape2;
		}
		else if(sym()==2){
			T nextshape2(locsusvs[i][0],locsusvs[i][1],fullsizes[i]);
			nextshape = nextshape2;
		}
		else{
			T nextshape2(locsusvs[i][0],locsusvs[i][1],
					locsusvs[i][2],fullsizes[i]);
			nextshape=nextshape2;
		}
		shapes.push_back(nextshape);
	}
	for(int i=0; i<N; i++) partvol += shapes[i].volume();
}

template<class T>
Torus<T>::Torus(Torus<T>& copy, int ndef){
	set_N(copy.get_N()); set_L(copy.get_L()); set_dim(copy.get_dim());
	partvol = copy.get_partvol(); ssizes = copy.get_sizes(); trng = copy.get_trng();
	std::vector<std::vector<gsl_vector*> > fs = ssizes.get_fullsizes();
	for(int i=0; i<N; i++){
		gsl_vector * pos = gsl_vector_alloc(dim);
		gsl_vector * u = gsl_vector_alloc(dim);
		gsl_vector * v = gsl_vector_alloc(dim);
		gsl_vector_memcpy(pos,copy.get_1_pos(i));
		if(sym()>1) gsl_vector_memcpy(u,copy.get_u(i));
		if(sym()>2) gsl_vector_memcpy(v,copy.get_v(i));
		T nextshape;
		if(sym()==1){T nextshape2(pos,fs[i]); nextshape = nextshape2;}
		else if(sym()==2){
			T nextshape2(pos, u, fs[i]);
			nextshape = nextshape2;
		}
		else{//sym==3 in this case
			T nextshape2(pos,u,v,fs[i]);
			nextshape = nextshape2;
		}
		shapes.push_back(nextshape);
	}
}

template<class T>
void Torus<T>::populate() {
	std::vector<std::vector<gsl_vector*> > fs = ssizes.get_fullsizes();
	for(int i=0; i<N; i++){
		gsl_vector * pos = gsl_vector_alloc(dim);
		gsl_vector * u = gsl_vector_alloc(dim);
		gsl_vector * v = gsl_vector_alloc(dim);
		for(int j=0; j<dim; j++){
			gsl_vector_set(pos,j,trng.num());
			if(sym()>1) gsl_vector_set(u,j,trng.num());
			if(sym()>2) gsl_vector_set(v,j,trng.num());
		}
		T nextshape;
		if(sym()==1){T nextshape2(pos,fs[i]); nextshape = nextshape2;}
		else if(sym()==2){
			T nextshape2(pos, u, fs[i]);
			nextshape = nextshape2;
			nextshape.normalize();
		}
		else{//sym==3 in this case
			T nextshape2(pos,u,v,fs[i]);
			nextshape = nextshape2;
			nextshape.normalize();
		}
		shapes.push_back(nextshape);
	}
}

template<class T>
void Torus<T>::find_set_L(double phi) {
	double V = partvol/phi;
	L = pow(V,1.0/(double)dim);
}

template<class T>
void Torus<T>::set_N(int N) {
	this->N = N;
}

template<class T>
void Torus<T>::set_L(double L) {
	this->L = L;
}

template<class T>
void Torus<T>::set_dim(int dim){
	this->dim = dim;
}

template<class T>
void Torus<T>::set_pos_coord(int i, int j, double x) {
	shapes[i].set_pos_coord(j,x);
}

template<class T>
void Torus<T>::set_u_coord(int i, int j, double x) {
	shapes[i].set_u_coord(j,x);
}

template<class T>
void Torus<T>::set_v_coord(int i, int j, double x) {
	shapes[i].set_v_coord(j,x);
}

template<class T>
double Torus<T>::R(int i, int j){
	return get_r(i)+get_r(j);
}

template<class T>
double Torus<T>::I(int i){
	return shapes[i].I();
}

template<class T>
gsl_vector * Torus<T>::F_loc(int i, int j, int k){
	return shapes[i].F_loc(shapes[j],k,L);
}

template<class T>
gsl_vector * Torus<T>::ell_1_vec(int i, int j, int k){
	return shapes[i].ell_vec(shapes[j],k,L);
}

template<class T>
gsl_vector * Torus<T>::ell_vec(int i, int j, int k){
	gsl_vector * reld = shapes[i].ell_vec(shapes[j],k,L);
	for(int i=0; i<dim; i++){
		gsl_vector_set(reld,i,L*gsl_vector_get(reld,i));
	}
	return reld;
}

template<class T>
double Torus<T>::ell2(int i, int j, int k) {
	return L*L*shapes[i].ell2(shapes[j], k,L);
}

template<class T>
double Torus<T>::pack_frac() {
	return partvol/std::pow(L,dim);
}

template<class T>
inline void Torus<T>::normalize(int i) {
	shapes[i].normalize();
}

template<class T>
inline void Torus<T>::modify(int i){
	for(int j=0; j<dim; j++){
		shapes[i].set_pos_coord(j,shapes[i].get_pos_coord(j));
	}
}

#endif /* TORUS_H_ */
