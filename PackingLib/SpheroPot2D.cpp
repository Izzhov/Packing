/*
 * SpheroPot2D.cpp
 *
 *  Created on: May 29, 2015
 *      Author: izzhov
 */

#include "SpheroPot2D.h"

SpheroPot2D::SpheroPot2D(SpheroTorus2D& ptorus):torus(ptorus) {
	epstrack=0;
	eps = 0.1;
	d = 0.000001; //for is_done
	kk = 1; P = 0.0001; //default values
	U=0; DL = 0; std::vector<double> DUi(2,0);
	std::vector<std::vector<double> > DUt(torus.get_N(),DUi); DU = DUt;
	std::vector<double> Dthi(torus.get_N(),0); Dtau = Dthi;
	calc_U_F();
}

void SpheroPot2D::next(){
	move(); calc_U_F();
}

void SpheroPot2D::calc_U_F() { //note forces are negatived *here*
	double L = torus.get_L();
	U = P*L*L; DL = -2*P*torus.get_L();//2PL wants smaller
	for(int i=0; i<torus.get_N(); i++){
		DU.at(i).at(0)=0; DU.at(i).at(1)=0; Dtau.at(i) = 0;
		for(int j=0; j<torus.get_N(); j++)
			for(int k=1; k<=4; k++) calc_U_F(i,j,k);
	}
}

void SpheroPot2D::calc_U_F(int i, int j, int k){
	if(i==j) return;
	double Rij = torus.R(i,j); double lijk2 = torus.ell2(i,j,k);
	if(lijk2>(Rij*Rij)) return;
	double lijk = sqrt(lijk2);
	if(i<j){//store energy
		U+=0.5*kk*(1-(lijk/Rij))*(1-(lijk/Rij));
	}
	double dudl = dUdl(lijk, Rij); double L = torus.get_L();
	double lis = torus.lis(i,j,k); double ljs = torus.ljs(i,j,k);
	double lijkx = L*torus.rel_x(i,j,k)+ljs*torus.get_u(j)-lis*torus.get_u(i);
	double lijky = L*torus.rel_y(i,j,k)+ljs*torus.get_v(j)-lis*torus.get_v(i);
	DU.at(i).at(0)+=-dudl*lijkx/lijk; //forces
	DU.at(i).at(1)+=-dudl*lijky/lijk;
	Dtau.at(i) += -(lis/lijk)*dudl*(torus.get_u(i)*lijky-torus.get_v(i)*lijkx);
	DL+=dudl/L; //wants bigger
}

double SpheroPot2D::dUdl(double lijk, double Rij){
	return (kk/Rij)*(1-(lijk/Rij))*lijk;
}

void SpheroPot2D::move() {
	translate(); rotate();
}

void SpheroPot2D::translate() {
	double L = torus.get_L();
	for(int i=0; i<torus.get_N(); i++){
		torus.set_1x(i,torus.get_1x(i)+(eps/L)*DU.at(i).at(0));
		torus.set_1y(i,torus.get_1y(i)+(eps/L)*DU.at(i).at(1));
	}
	torus.set_L(L+eps*DL);//2*P*L wants smaller
}

void SpheroPot2D::rotate() { //needs fixin for theta!
	double uu; double vv; double theta;
	for(int i=0; i<torus.get_N();i++){
		uu = torus.get_u(i); vv = torus.get_v(i);
		theta = (eps/torus.Iom(i))*Dtau[i];
		torus.set_u(i,cos(theta)*uu-sin(theta)*vv);
		torus.set_v(i,sin(theta)*uu+cos(theta)*vv);
		torus.normalize(i);
	}

}

double SpheroPot2D::next_eps() {
	move();
	epstrack=(calc_U()-U)/D_mag_2();
	calc_U_F();
	return epstrack;
}

double SpheroPot2D::D_mag_2() {
	int N = torus.get_N();
	double mag = 0; std::vector<double> q; double x; double y;
	for(std::vector<std::vector<double> >::iterator it = DU.begin();
			it != DU.end(); ++it){
		q = *it; x = q.at(0); y = q.at(1); mag-=x*x+y*y;
	}
	for(int i = 0; i<N; i++){
		mag-=(Dtau[i]*Dtau[i])/torus.Iom(i);
	}
	mag-=DL*DL;
	return mag;
}

double SpheroPot2D::calc_U() {
	double L = torus.get_L(); double newU=P*L*L;
	for(int i=0; i<torus.get_N(); i++)
		for(int j=0; j<torus.get_N(); j++)
			for(int k=1; k<=4; k++) newU+=calc_U(i,j,k);
	return newU;
}

double SpheroPot2D::calc_U(int i, int j, int k) {
	double calcU=0;
	if(i==j) return 0;
	double Rij = torus.R(i,j); double lijk2 = torus.ell2(i,j,k);
	if(lijk2>(Rij*Rij)) return 0;
	double lijk = sqrt(lijk2);
	if(i<j){//store energy
		calcU+=0.5*kk*(1-(lijk/Rij))*(1-(lijk/Rij));
	}
	return calcU;
}
