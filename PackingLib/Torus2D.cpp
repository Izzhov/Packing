/*
 * Torus2D.cpp
 *
 *  Created on: May 21, 2015
 *      Author: root
 */

#include "Torus2D.h"

Torus2D::Torus2D(){
	set_seed(trng.get_seed());
	N = 0; L = 1; partarea = 0;
}

Torus2D::Torus2D(int N, std::vector<double> sizes, std::vector<double> wts) {
	Torus2D intmd(N,1,sizes,wts);
	*this = intmd;
}

Torus2D::Torus2D(int N, std::vector<double> sizes, std::vector<double> wts, int t){
	TimeRNG01 newrng(t); trng=newrng;
	set_seed(t);
	partarea=0;
	set_N(N);
	set_L(1);
	ssizes.construct(N,sizes,wts);
	populate();
	for(int i = 0; i<N; i++) partarea += sphs.at(i).area();
}

Torus2D::Torus2D(int N, double L, std::vector<double> sizes, std::vector<double> wts) {
	set_seed(trng.get_seed());
	partarea = 0;
	set_N(N);
	set_L(L);
	ssizes.construct(N, sizes, wts);
	populate();
	for(int i = 0; i<N; i++) partarea += sphs.at(i).area();
}

Torus2D::Torus2D(int N, double L, std::vector<std::vector<double> > locs){
	set_seed(0);
	partarea=0;
	set_N(N); set_L(L);
	set_populate(locs);
	for(int i=0; i<N; i++) partarea+=sphs.at(i).area();
}

void Torus2D::set_seed(int seed){
	this->seed = seed;
}

void Torus2D::populate() {
	double x; double y; double r;
	for(int i = 0; i<N; i++){
		x = trng.num(); y = trng.num();
		r = ssizes.get_radii().at(i);
		Sphere2D nextsph(x,y,r);
		sphs.push_back(nextsph);
	}
}

void Torus2D::set_populate(std::vector<std::vector<double> > locs){
	double x; double y;
	for(int i = 0; i<N; i++){
		x = locs.at(i).at(0); y = locs.at(i).at(1);
		Sphere2D nextsph(x,y,1); sphs.push_back(nextsph);
	}
}

void Torus2D::find_set_L(double phi){
	double A = partarea/phi;
	L = sqrt(A);
}

void Torus2D::set_N(int N) {
	this->N = N;
}

void Torus2D::set_L(double L) {
	this->L = L;
}

void Torus2D::set_1x(int i, double x){sphs.at(i).set_x(x);}
void Torus2D::set_1y(int i, double y){sphs.at(i).set_y(y);}

double Torus2D::ell(int i, int j, int k) {
	double rx = rel_x(i,j,k); double ry = rel_y(i,j,k);
	return sqrt(rx*rx+ry*ry);
}

double Torus2D::pack_frac() {
	return partarea/(L*L);
}
