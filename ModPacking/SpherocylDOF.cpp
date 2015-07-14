/*
 * SpherocylDOF.cpp
 *
 *  Created on: Jul 6, 2015
 *      Author: izzhov
 */

#include "SpherocylDOF.h"

SpherocylDOF::SpherocylDOF() {}

SpherocylDOF::SpherocylDOF(gsl_vector* pos, std::vector<gsl_vector*> size) {}

SpherocylDOF::SpherocylDOF(gsl_vector* pos, gsl_vector* u,
		std::vector<gsl_vector*> size) {
	set_pos(pos); set_u(u); set_r(gsl_vector_get(size[0],0));
	set_a(gsl_vector_get(size[1],0));
}

SpherocylDOF::SpherocylDOF(gsl_vector* pos, gsl_vector* u, gsl_vector* v,
		std::vector<gsl_vector*> size) {}

SpherocylDOF::SpherocylDOF(SpherocylDOF * orig, double ratio){
	int size = (orig->pos)->size;
	gsl_vector * newpos = gsl_vector_alloc(size);
	gsl_vector * newu = gsl_vector_alloc(size);
	gsl_vector_memcpy(newpos,orig->pos);
	gsl_vector_memcpy(newu,orig->u);
	set_pos(newpos); set_u(newu);
	set_r(ratio*(*orig).get_r());
	set_a((1/ratio)*((*orig).get_a()-1)+1);
}

void SpherocylDOF::set_pos(gsl_vector* pos) {
	this->pos = pos;
}

void SpherocylDOF::set_u(gsl_vector* u) {
	this->u = u;
}

void SpherocylDOF::set_r(double r) {
	this->r = r;
}

void SpherocylDOF::set_a(double a) {
	this->a = a;
}

void SpherocylDOF::set_pos_coord(int i, double x) {
	gsl_vector_set(pos,i,mm::pmod(x,1));
}

void SpherocylDOF::set_u_coord(int i, double x) {
	gsl_vector_set(u,i,x);
}

double SpherocylDOF::sphere_volume(){
	if(pos->size == 2) return (M_PI*r*r);
	else return (4.0/3.0)*M_PI*r*r*r;
}

double SpherocylDOF::cyl_volume(){
	if(pos->size == 2) return 4*(a-1)*r*r;
	else return M_PI*r*r*(2*(a-1)*r);
}

double SpherocylDOF::volume() {
	return sphere_volume()+cyl_volume();
}

double SpherocylDOF::I() {
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

double SpherocylDOF::max_d() {
	return (a-1)*r;
}

double SpherocylDOF::lsl(){
	return 2*a*r;
}

gsl_vector* SpherocylDOF::lisljs(SpherocylDOF s, int k, double L, int ncon) {
	gsl_vector * thels = gsl_vector_calloc(2);
	double lis = 0; double ljs = 0;
	double li = 2*max_d(); double lj = 2*s.max_d();
	double uiuj; gsl_blas_ddot(get_u(),s.get_u(),&uiuj);
	gsl_vector * rij = mm::rel(pos,s.get_pos(),k);
	gsl_vector_scale(rij,L);
	double uirij; gsl_blas_ddot(get_u(), rij, &uirij);
	double ujrij; gsl_blas_ddot(s.get_u(),rij,&ujrij);
	int si = mm::sgn(uirij); int sj = mm::sgn(ujrij);
	if((1-uiuj*uiuj)<std::pow(10,-8)){//then spherocyls are parallel
		if(std::abs(uirij)>(li+lj)/2.0){//force only applied in one place
			if(ncon==0){
				lis = si*li/2.0; ljs = -sj*lj/2.0;
			}
			//else pick the ones that'll never be touchin'
			else{
				lis = -si*li/2.0; ljs = sj*lj/2.0;
			}
		}
		else{//force applied at 2 forcends (see6-16-15 page 2 for algo)
			if(ncon==0){
				lis = mm::minmag(si,si*li/2.0,uirij+si*lj/2.0);
				ljs = mm::minmag(sj,sj*lj/2.0,-ujrij+sj*li/2.0);
			}
			else if(ncon==1){
				lis = mm::minmag(-si,-si*li/2.0,uirij-si*lj/2.0);
				ljs = mm::minmag(-sj,-sj*lj/2.0,-ujrij-sj*li/2.0);
			}
			else{
				double posi = std::min(li/2.0,uirij+lj/2.0);
				double negi = std::max(-li/2.0,uirij-lj/2.0);
				double posj = std::min(lj/2.0,-ujrij+li/2.0);
				double negj = std::max(-lj/2.0,-ujrij-li/2.0);
				lis = (posi+negi)/2.0; ljs = (posj+negj)/2.0;
			}
		}
	}
	else{//not parallel, use regular algorithm from paper
		double lip = (uirij-uiuj*ujrij)/(1-uiuj*uiuj);
		double ljp = (uiuj*uirij-ujrij)/(1-uiuj*uiuj);
		double bli = std::abs(lip)-li/2.0;
		double blj = std::abs(ljp)-lj/2.0;
		double bbli = std::abs(lip)+li/2.0;
		double bblj = std::abs(ljp)+lj/2.0;
		if(ncon==0){
			if(blj>=bli){
				ljs = mm::sign(lj/2.0,ljp);
				lis = std::max(-li/2.0,std::min(uirij+ljs*uiuj,li/2.0));
			}
			else{
				lis = mm::sign(li/2.0,lip);
				ljs = std::max(-lj/2.0,std::min(-ujrij+lis*uiuj,lj/2.0));
			}
		}
		else if(ncon==1){
			double lim=0; double ljm=0;
			if(blj>=bli){
				ljm = mm::sign(lj/2.0,ljp);
				lim = std::max(-li/2.0,std::min(uirij+ljm*uiuj,li/2.0));
			}
			else{
				lim = mm::sign(li/2.0,lip);
				ljm = std::max(-lj/2.0,std::min(-ujrij+lim*uiuj,lj/2.0));
			}
			if(bblj<=bbli){
				ljs=mm::sign(lj/2.0,-ljp);
				lis=std::max(-li/2.0,std::min(uirij+ljs*uiuj,li/2.0));
			}
			else{
				lis=mm::sign(li/2.0,-lip);
				ljs=std::max(-lj/2.0,std::min(-ujrij+lis*uiuj,lj/2.0));
			}
			if(std::abs(lis-lim)<std::pow(10,-8) || std::abs(ljs-ljm)<std::pow(10,-8)){
				lis=0; ljs=std::pow(10,5);
			}
		}
		else if(ncon==2){
			if(bli<=0 && blj <=0){lis = lip; ljs = ljp;}
			else{lis=0; ljs=std::pow(10,5);}
		}
	}
	gsl_vector_set(thels,0,lis); gsl_vector_set(thels,1,ljs);
	gsl_vector_free(rij);
	return thels;
}

gsl_vector* SpherocylDOF::F_loc(SpherocylDOF s, int k, double L, int ncon) {
	gsl_vector * floc = gsl_vector_alloc(pos->size);
	gsl_vector * thels = lisljs(s,k,L,ncon);
	double li = gsl_vector_get(thels,0);
	gsl_vector_memcpy(floc,get_u()); gsl_vector_scale(floc,li);
	gsl_vector_free(thels);
	return floc;
}

gsl_vector* SpherocylDOF::ell_vec(SpherocylDOF s, int k, double L, int ncon) {
	gsl_vector * rij = mm::rel(pos,s.get_pos(),k);
	gsl_vector * uu = gsl_vector_alloc(rij->size);
	gsl_vector * thels = lisljs(s,k,L,ncon);
	double li = gsl_vector_get(thels,0); double lj = gsl_vector_get(thels,1);
	gsl_vector_memcpy(uu,get_u()); gsl_vector_scale(uu,-li/L);
	gsl_vector_add(rij,uu);
	gsl_vector_memcpy(uu,s.get_u()); gsl_vector_scale(uu,lj/L);
	gsl_vector_add(rij,uu);
	gsl_vector_free(uu); gsl_vector_free(thels);
	return rij;
}

double SpherocylDOF::ell2(SpherocylDOF s, int k, double L, int ncon) {
	double dd;//the distance
	gsl_vector * reld = ell_vec(s,k,L, ncon);
	gsl_blas_ddot(reld,reld,&dd);
	gsl_vector_free(reld);
	return dd;
}

gsl_vector* SpherocylDOF::lisljs(SpherocylDOF s, double L, int ncon) {
	gsl_vector * thels = gsl_vector_calloc(2);
	double lis = 0; double ljs = 0;
	double li = 2*max_d(); double lj = 2*s.max_d();
	double uiuj; gsl_blas_ddot(get_u(),s.get_u(),&uiuj);
	gsl_vector * rij = mm::rel(pos,s.get_pos());
	gsl_vector_scale(rij,L);
	double uirij; gsl_blas_ddot(get_u(), rij, &uirij);
	double ujrij; gsl_blas_ddot(s.get_u(),rij,&ujrij);
	int si = mm::sgn(uirij); int sj = mm::sgn(ujrij);
	if((1-uiuj*uiuj)<std::pow(10,-8)){//then spherocyls are parallel
		if(std::abs(uirij)>(li+lj)/2.0){//force only applied in one place
			if(ncon==0){
				lis = si*li/2.0; ljs = -sj*lj/2.0;
			}
			//else pick the ones that'll never be touchin'
			else{
				lis = -si*li/2.0; ljs = sj*lj/2.0;
			}
		}
		else{//force applied at 2 forcends (see6-16-15 page 2 for algo)
			if(ncon==0){
				lis = mm::minmag(si,si*li/2.0,uirij+si*lj/2.0);
				ljs = mm::minmag(sj,sj*lj/2.0,-ujrij+sj*li/2.0);
			}
			else if(ncon==1){
				lis = mm::minmag(-si,-si*li/2.0,uirij-si*lj/2.0);
				ljs = mm::minmag(-sj,-sj*lj/2.0,-ujrij-sj*li/2.0);
			}
			else{
				double posi = std::min(li/2.0,uirij+lj/2.0);
				double negi = std::max(-li/2.0,uirij-lj/2.0);
				double posj = std::min(lj/2.0,-ujrij+li/2.0);
				double negj = std::max(-lj/2.0,-ujrij-li/2.0);
				lis = (posi+negi)/2.0; ljs = (posj+negj)/2.0;
			}
		}
	}
	else{//not parallel, use regular algorithm from paper
		double lip = (uirij-uiuj*ujrij)/(1-uiuj*uiuj);
		double ljp = (uiuj*uirij-ujrij)/(1-uiuj*uiuj);
		double bli = std::abs(lip)-li/2.0;
		double blj = std::abs(ljp)-lj/2.0;
		double bbli = std::abs(lip)+li/2.0;
		double bblj = std::abs(ljp)+lj/2.0;
		if(ncon==0){
			if(blj>=bli){
				ljs = mm::sign(lj/2.0,ljp);
				lis = std::max(-li/2.0,std::min(uirij+ljs*uiuj,li/2.0));
			}
			else{
				lis = mm::sign(li/2.0,lip);
				ljs = std::max(-lj/2.0,std::min(-ujrij+lis*uiuj,lj/2.0));
			}
		}
		else if(ncon==1){
			double lim=0; double ljm=0;
			if(blj>=bli){
				ljm = mm::sign(lj/2.0,ljp);
				lim = std::max(-li/2.0,std::min(uirij+ljm*uiuj,li/2.0));
			}
			else{
				lim = mm::sign(li/2.0,lip);
				ljm = std::max(-lj/2.0,std::min(-ujrij+lim*uiuj,lj/2.0));
			}
			if(bblj<=bbli){
				ljs=mm::sign(lj/2.0,-ljp);
				lis=std::max(-li/2.0,std::min(uirij+ljs*uiuj,li/2.0));
			}
			else{
				lis=mm::sign(li/2.0,-lip);
				ljs=std::max(-lj/2.0,std::min(-ujrij+lis*uiuj,lj/2.0));
			}
			if(std::abs(lis-lim)<std::pow(10,-8) || std::abs(ljs-ljm)<std::pow(10,-8)){
				lis=0; ljs=std::pow(10,5);
			}
		}
		else if(ncon==2){
			if(bli<=0 && blj <=0){lis = lip; ljs = ljp;}
			else{lis=0; ljs=std::pow(10,5);}
		}
	}
	gsl_vector_set(thels,0,lis); gsl_vector_set(thels,1,ljs);
	gsl_vector_free(rij);
	return thels;
}

gsl_vector* SpherocylDOF::F_loc(SpherocylDOF s, double L, int ncon) {
	gsl_vector * floc = gsl_vector_alloc(pos->size);
	gsl_vector * thels = lisljs(s,L,ncon);
	double li = gsl_vector_get(thels,0);
	gsl_vector_memcpy(floc,get_u()); gsl_vector_scale(floc,li);
	gsl_vector_free(thels);
	return floc;
}

gsl_vector* SpherocylDOF::ell_vec(SpherocylDOF s, double L, int ncon) {
	gsl_vector * rij = mm::rel(pos,s.get_pos());
	gsl_vector * uu = gsl_vector_alloc(rij->size);
	gsl_vector * thels = lisljs(s,L,ncon);
	double li = gsl_vector_get(thels,0); double lj = gsl_vector_get(thels,1);
	gsl_vector_memcpy(uu,get_u()); gsl_vector_scale(uu,-li/L);
	gsl_vector_add(rij,uu);
	gsl_vector_memcpy(uu,s.get_u()); gsl_vector_scale(uu,lj/L);
	gsl_vector_add(rij,uu);
	gsl_vector_free(uu); gsl_vector_free(thels);
	return rij;
}

double SpherocylDOF::ell2(SpherocylDOF s, double L, int ncon) {
	double dd;//the distance
	gsl_vector * reld = ell_vec(s,L, ncon);
	gsl_blas_ddot(reld,reld,&dd);
	gsl_vector_free(reld);
	return dd;
}

bool SpherocylDOF::touch(SpherocylDOF s, double L){
	bool dotheytouch = false;
	for(int k=1; k<=mm::int_pow(2,pos->size);k++){
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
		double ll2; gsl_blas_ddot(rij,rij,&ll2);
		ll2+=2*ujrij*ljs-2*uirij*lis-2*uiuj*lis*ljs+ljs*ljs+lis*lis;
		if(std::sqrt(ll2)<=get_r()+s.get_r()){
			dotheytouch = true;
		}
		gsl_vector_free(rij);
	}
	return dotheytouch;
}

bool SpherocylDOF::is_far(SpherocylDOF s, double L){
	bool isfar = false;
	double lis = 0; double ljs = 0;
	double li = 2*max_d(); double lj = 2*s.max_d();
	double uiuj; gsl_blas_ddot(get_u(),s.get_u(),&uiuj);
	gsl_vector * rij = mm::rel(pos,s.get_pos());
	gsl_vector_scale(rij,L);
	double uirij; gsl_blas_ddot(get_u(), rij, &uirij);
	double ujrij; gsl_blas_ddot(s.get_u(),rij,&ujrij);
	lis=li/2.0; ljs=lj/2.0;
	double ll1; gsl_blas_ddot(rij,rij,&ll1);
	ll1+=2*ujrij*ljs-2*uirij*lis-2*uiuj*lis*ljs+ljs*ljs+lis*lis;
	ll1=std::sqrt(ll1);
	lis=-li/2.0; ljs=-lj/2.0;
	double ll2; gsl_blas_ddot(rij,rij,&ll2);
	ll2+=2*ujrij*ljs-2*uirij*lis-2*uiuj*lis*ljs+ljs*ljs+lis*lis;
	ll2=std::sqrt(ll2);
	isfar=(ll1>=s.get_r()-get_r())||(ll2>=s.get_r()-get_r());
	gsl_vector_free(rij);
	return isfar;
}

void SpherocylDOF::normalize() {
	double unorm = gsl_blas_dnrm2(u);
	gsl_vector_scale(u,1.0/unorm);
}

gsl_vector* SpherocylDOF::get_v() {}
