/*
 * SpheroTorus2D.cpp
 *
 *  Created on: May 27, 2015
 *      Author: izzhov
 */

#include "SpheroTorus2D.h"

SpheroTorus2D::SpheroTorus2D() {
	N=0; L=1; partarea=0;
}

SpheroTorus2D::SpheroTorus2D(int N, std::vector<double> rs,
		std::vector<double> as, std::vector<double> wts) {
	SpheroTorus2D intmd(N,1,rs,as,wts);
	*this = intmd;
}

SpheroTorus2D::SpheroTorus2D(int N, double L, std::vector<double> rs,
		std::vector<double> as, std::vector<double> wts) {
	partarea=0;
	set_N(N); set_L(L);
	ssizes.construct(N,rs,as,wts);
	populate();
	for(int i = 0; i<N; i++) partarea += sphs.at(i).area();
}

SpheroTorus2D::SpheroTorus2D(int N, double L,
		std::vector<std::vector<double> > locs) {
	partarea=0;
	set_N(N); set_L(L);
	set_populate(locs);
	for(int i=0; i<N; i++) partarea+=sphs.at(i).area();
}

void SpheroTorus2D::populate() {
	double x; double y; double th; double u; double v; double r; double a;
	for(int i = 0; i<N; i++){
		x = trng.num(); y = trng.num(); th = M_PI*trng.num();
		u = cos(th); v = sin(th);
		r = ssizes.get_radii().at(i); a = ssizes.get_lengths().at(i);
		Sphero2D nextsph(x,y,u,v,r,a); nextsph.normalize();
		sphs.push_back(nextsph);
	}
}

void SpheroTorus2D::set_populate(std::vector<std::vector<double> > locs) {
	double x; double y; double u; double v;
	for(int i = 0; i<N; i++){
		x = locs.at(i).at(0); y = locs.at(i).at(1);
		u = locs.at(i).at(2); v = locs.at(i).at(3);
		Sphero2D nextsph(x,y,u,v,1,2); nextsph.normalize();
		sphs.push_back(nextsph);
	}
}

void SpheroTorus2D::set_N(int N) {this->N = N;}

void SpheroTorus2D::set_L(double L) {
	this->L = L;
}

void SpheroTorus2D::set_1x(int i, double x) {sphs.at(i).set_x(x);}

void SpheroTorus2D::set_1y(int i, double y) {sphs.at(i).set_y(y);}

void SpheroTorus2D::set_u(int i, double u) {
	sphs.at(i).set_u(u);
}

void SpheroTorus2D::set_v(int i, double v) {
	sphs.at(i).set_v(v);
}

void SpheroTorus2D::find_set_L(double phi) {
	double A = partarea/phi;
	L = sqrt(A);
}

void SpheroTorus2D::normalize(int i){
	sphs.at(i).normalize();
}

double SpheroTorus2D::rel_x(int i, int j, int k) {
	if(k==1||k==4)
		return mymath::pmod(sphs.at(j).get_x()-sphs.at(i).get_x(),1);
	else
		return mymath::pmod(sphs.at(j).get_x()-sphs.at(i).get_x(),1)-1;
}

double SpheroTorus2D::rel_y(int i, int j, int k) {
	if(k==1||k==2)
		return mymath::pmod(sphs.at(j).get_y()-sphs.at(i).get_y(),1);
	else
		return mymath::pmod(sphs.at(j).get_y()-sphs.at(i).get_y(),1)-1;
}

double SpheroTorus2D::len(int i) {
	return 2*(sphs.at(i).get_a()-1)*sphs.at(i).get_r();
}

double SpheroTorus2D::near(int i, int j){
	return sphs.at(i).get_a()*sphs.at(i).get_r()+
			sphs.at(j).get_a()*sphs.at(j).get_r();
}

double SpheroTorus2D::R(int i, int j){
	return sphs.at(i).get_r()+sphs.at(j).get_r();
}

double SpheroTorus2D::Iom(int i){
	return sphs.at(i).Iom();
}

double SpheroTorus2D::rij2(int i, int j, int k) {
	double xx = rel_x(i,j,k); double yy = rel_y(i,j,k);
	return L*L*(xx*xx+yy*yy);
}

double SpheroTorus2D::ui_rij(int i, int j, int k) {
	double xx = rel_x(i,j,k); double yy = rel_y(i,j,k);
	double ui = sphs.at(i).get_u(); double vi = sphs.at(i).get_v();
	return L*(ui*xx+vi*yy);
}

double SpheroTorus2D::uj_rij(int i, int j, int k) {
	double xx = rel_x(i,j,k); double yy = rel_y(i,j,k);
	double uj = sphs.at(j).get_u(); double vj = sphs.at(j).get_v();
	return L*(uj*xx+vj*yy);
}

double SpheroTorus2D::ui_uj(int i, int j, int k) {
	double ui = sphs.at(i).get_u(); double vi = sphs.at(i).get_v();
	double uj = sphs.at(j).get_u(); double vj = sphs.at(j).get_v();
	return (ui*uj+vi*vj);
}

double SpheroTorus2D::lip(int i, int j, int k) {
	double ui_rijj = ui_rij(i,j,k); double uj_rijj = uj_rij(i,j,k);
	double ui_ujj = ui_uj(i,j,k);
	return (ui_rijj-(ui_ujj)*(uj_rijj))/(1-ui_ujj*ui_ujj);
}

double SpheroTorus2D::ljp(int i, int j, int k) {
	double ui_rijj = ui_rij(i,j,k); double uj_rijj = uj_rij(i,j,k);
	double ui_ujj = ui_uj(i,j,k);
	return (ui_ujj*ui_rijj-uj_rijj)/(1-ui_ujj*ui_ujj);
}

double SpheroTorus2D::bli(int i, int j, int k) {
	return std::abs(lip(i,j,k))-len(i)/2;
}

double SpheroTorus2D::blj(int i, int j, int k) {
	return std::abs(ljp(i,j,k))-len(j)/2;
}

double SpheroTorus2D::lis(int i, int j, int k) {
	double blii = bli(i,j,k); double bljj = blj(i,j,k);
	double lipp = lip(i,j,k);
	double li = len(i);
	double ui_rijj = ui_rij(i,j,k); double ui_ujj = ui_uj(i,j,k);
	if (blii<=0 && bljj<=0) return lipp;
	else if (bljj>=blii)
		return std::max(-li/2,std::min(ui_rijj+ljs(i,j,k)*ui_ujj,li/2));
	else return mymath::sign(li/2,lipp);
}

double SpheroTorus2D::ljs(int i, int j, int k) {
	double blii = bli(i,j,k); double bljj = blj(i,j,k);
	double ljpp = ljp(i,j,k);
	double lj = len(j);
	double uj_rijj = uj_rij(i,j,k); double ui_ujj = ui_uj(i,j,k);
	if (blii <=0 && bljj<=0) return ljpp;
	else if (bljj>=blii) return mymath::sign(lj/2, ljpp);
	else return std::max(-lj/2,std::min(-uj_rijj+lis(i,j,k)*ui_ujj,lj/2));
}

double SpheroTorus2D::ell2(int i, int j, int k) {
	double rij22 = rij2(i,j,k);
	if(rij22 > near(i,j)*near(i,j)) return rij22;
	//^not accurate, but it doesn't matter for now
	double ui_ujj = ui_uj(i,j,k);
	double ui_rijj = ui_rij(i,j,k);
	if ((1-(ui_ujj*ui_ujj))<std::pow(10,-8)){
		double li = len(i); double lj = len(j);
		double maxx = std::max(0.0,std::abs(ui_rijj)-(li+lj)/2);
		return rij22-ui_rijj*ui_rijj+maxx*maxx;
	}
	else{
		double liss = lis(i,j,k); double ljss = ljs(i,j,k);
		double uj_rijj = uj_rij(i,j,k);
		return rij22 + 2*uj_rijj*ljss-2*ui_rijj*liss-2*ui_ujj*liss*ljss
				+ljss*ljss+liss*liss;
	}

}

double SpheroTorus2D::pack_frac() {
	return partarea/(L*L);
}
