/*
 * Potential2DNoNbr.cpp
 *
 *  Created on: May 26, 2015
 *      Author: izzhov
 */

#include "Potential2DNoNbr.h"

//constructor for harm. osc's
Potential2DNoNbr::Potential2DNoNbr(Torus2D& ptorus):torus(ptorus){
	eps = 0.1; //currently: L^1 for dividin translation
	d = 0.000001; //for is_done
	kk = 1; P = std::pow(10,-4); //default values
	U=0; DL = 0; std::vector<double> DUi(2,0);
	std::vector<std::vector<double> > DUt(torus.get_N(),DUi); DU = DUt;
	harm_osc(); grad_harm_osc(); L_harm_osc();
}

double Potential2DNoNbr::D_mag_2(){
	double mag = 0; std::vector<double> q; double x; double y;
	for(std::vector<std::vector<double> >::iterator it = DU.begin();
			it != DU.end(); ++it){
		q = *it; x = q.at(0); y = q.at(1); mag+=x*x+y*y;
	}
	mag+=DL*DL;
	return mag;
}

void Potential2DNoNbr::harm_osc() {
	U = 0; double L = torus.get_L(); int N = torus.get_N();
	U += P*L*L; double lijk; double Rij;
	for(int i = 0; i < N; i++){
		for(int j=i+1; j < N; j++){
			for(int k = 1; k <= 4; k++){
				lijk = torus.ell(i,j,k); Rij = torus.R(i,j);
				if(lijk<=Rij){
					U += 0.5*kk*(1-(lijk/Rij))*(1-(lijk/Rij));
				}
			}
		}
	}
}

void Potential2DNoNbr::grad_harm_osc(){
	int N = torus.get_N(); double lijk; double Rij;
	for(int i=0; i<N; i++){
		DU.at(i).at(0)=0; DU.at(i).at(1)=0;
		for(int j=0; j<i; j++){ //switch i and j for this!
			Rij = torus.R(j,i);
			for(int k = 1; k<=4; k++){
				lijk = torus.ell(j,i,k);
				if(lijk<=Rij){
					DU.at(i).at(0)+=(kk/Rij)*((lijk/Rij)-1)*dxj_ell(j,i,k);
					DU.at(i).at(1)+=(kk/Rij)*((lijk/Rij)-1)*dyj_ell(j,i,k);
				}
			}
		}
		for(int j = i+1; j<N; j++){ //don't switch
			Rij = torus.R(i,j);
			for(int k=1; k<=4; k++){
				lijk = torus.ell(i,j,k);
				if(lijk<=Rij){
					DU.at(i).at(0)+=(kk/Rij)*((lijk/Rij)-1)*dxi_ell(i,j,k);
					DU.at(i).at(1)+=(kk/Rij)*((lijk/Rij)-1)*dyi_ell(i,j,k);
				}
			}
		}
	}
}

double Potential2DNoNbr::dxj_ell(int i, int j, int k){
	double L = torus.get_L();
	return L*torus.rel_x(i,j,k)/torus.ell(i,j,k);
}

double Potential2DNoNbr::dxi_ell(int i, int j, int k){
	return -dxj_ell(i,j,k);
}

double Potential2DNoNbr::dyj_ell(int i, int j, int k){
	double L = torus.get_L();
	return L*torus.rel_y(i,j,k)/torus.ell(i,j,k);
}

double Potential2DNoNbr::dyi_ell(int i, int j, int k){
	return -dyj_ell(i,j,k);
}

void Potential2DNoNbr::L_harm_osc(){
	DL = 0; int N = torus.get_N();
	DL += 2*P*torus.get_L(); double lijk; double Rij;
	for(int i = 0; i<N; i++){
		for(int j = i+1; j<N; j++){
			for(int k = 1; k<=4; k++){
				lijk = torus.ell(i,j,k); Rij = torus.R(i,j);
				if(lijk<=Rij) DL+=(kk/Rij)*((lijk/Rij)-1)*dL_ell(i,j,k);
			}
		}
	}
}

double Potential2DNoNbr::dL_ell(int i, int j, int k){
	return torus.ell(i,j,k)/torus.get_L();
}

void Potential2DNoNbr::full_next(){
	move(); harm_osc(); grad_harm_osc(); L_harm_osc();
}

void Potential2DNoNbr::full_next_debug(){
	next_epscheck(); next_rest();
}

double Potential2DNoNbr::next_epscheck(){
	double oldU = U; move(); harm_osc();
	return std::abs(U-oldU)/D_mag_2();
}

void Potential2DNoNbr::next_rest(){
	grad_harm_osc(); L_harm_osc();
}

void Potential2DNoNbr::move(){
	double L = torus.get_L(); //for incrementing x by smaller
	for(int i = 0; i<torus.get_N(); i++){
		torus.set_1x(i,torus.get_1x(i)-(eps/(L))*DU.at(i).at(0));
		torus.set_1y(i,torus.get_1y(i)-(eps/(L))*DU.at(i).at(1));
	}
	torus.set_L(L-eps*DL);
}

bool Potential2DNoNbr::is_done(){
	bool id = true; int N = torus.get_N(); int L = torus.get_L();
	double force; double fx; double fy;
	for(int i=0; i<N; i++){
		fx = DU.at(i).at(0); fy = DU.at(i).at(1);
		force = std::sqrt(fx*fx+fy*fy);
		if(force >= P*d){
			id = false; break;
		}
	}
	return (id && -2*P*L<DL && DL<P*L);
}
