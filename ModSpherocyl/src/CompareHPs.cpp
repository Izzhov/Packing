/*
 * CompareHPs.cpp
 *
 *  Created on: Jul 10, 2015
 *      Author: izzhov
 */

#include "CompareHPs.h"

CompareHPs::CompareHPs(HarmPot<Torus<SpherocylDOF> >& ppot,
		Torus<SpherocylDOF>& bbox,
		HarmPotNbrList<Torus<SpherocylDOF>, SpherocylDOF>& ppot2,
		Torus<SpherocylDOF>& bbox2):pot(ppot),box(bbox),pot2(ppot2),box2(bbox2){
	alph=0; alph2=alph; eps=std::pow(10,-4); eps2=eps; sig=0.1; imax = 1000000;
}

bool CompareHPs::minimize() {
	std::cout << std::setprecision(14);
	std::cout << "L1: " << box.get_L() << std::endl;
	std::cout << "L2: " << box2.get_L() << std::endl;
	int i=0; double del; double del0; double del2; double del02;
	double eta; double etaprev; double eta2; double etaprev2;
	//r stores negative deriv; rsec stores pos deriv after moving sig0
	gsl_vector * r = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	gsl_vector * rsec = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	gsl_vector * r2 = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	gsl_vector * rsec2 = gsl_vector_alloc(box.get_N()*box.get_dim()*box.sym()+1);
	//r=b-Ax
	pot.calc_U_F(); pot2.calc_U_F();
	std::cout << "DL1: " << pot.get_DL() << std::endl;
	std::cout << "DL2: " << pot2.get_DL() << std::endl;
	store_deriv(r,-1,1); store_deriv(r2,-1,2);
	//del=rTr
	del = unit_fix(r,1); del2 = unit_fix(r2,2);
	del0 = del; del02 = del2;
	if(del!=del2) std::cout << "messed up" << std::endl;
	while(i<imax && del>eps*eps*del0){
		alph=-sig; alph2=-sig;
		move(sig,1);
		//std::cout << "L1: " << box.get_L() << std::endl;
		move(sig,2);
		//std::cout << "L2: " << box2.get_L() << std::endl;
		pot.calc_U_F();
		//std::cout << "DL1: " << pot.get_DL() << std::endl;
		pot2.calc_U_F();
		//std::cout << "DL2: " << pot2.get_DL() << std::endl;
		store_deriv(rsec,1,1); store_deriv(rsec2,1,2);
		gsl_blas_ddot(rsec,r,&etaprev); gsl_blas_ddot(rsec2,r2,&etaprev2);
		gsl_blas_ddot(r,r,&eta); gsl_blas_ddot(r2,r2,&eta2);
		eta=-eta; eta2=-eta2;
		alph = alph*eta/(etaprev-eta);
		//std::cout << "alph: " << alph << std::endl;
		alph2 = alph2*eta2/(etaprev2-eta2);
		//std::cout << "alph2: " << alph2 << std::endl;
		if(alph>100) alph=100;
		if(alph>=sig){//only do second move if it's in the right direction
			//need to load old derivs before second move
			load_deriv(r,-1,1);
			move(alph-sig,1);
			//std::cout << "L1: " << box.get_L() << std::endl;
			pot.calc_U_F(); store_deriv(r,-1,1);
			//std::cout << "DL1: " << pot.get_DL() << std::endl;
		}
		if(alph2>100) alph2=100;
		if(alph2>=sig){//only do second move if it's in the right direction
			//need to load old derivs before second move
			load_deriv(r2,-1,2);
			move(alph2-sig,2);
			//std::cout << "L2: " << box2.get_L() << std::endl;
			pot2.calc_U_F(); store_deriv(r2,-1,2);
			//std::cout << "DL2: " << pot.get_DL() << std::endl;
		}
		del = unit_fix(r,1); del2 = unit_fix(r2,2);
		del0 = box.get_dim()*box.get_dim()*pot.get_P()*pot.get_P()*
				mm::int_pow(box.get_L(),2*box.get_dim()-2);
		del02 = box.get_dim()*box.get_dim()*pot.get_P()*pot.get_P()*
				mm::int_pow(box2.get_L(),2*box.get_dim()-2);
		if(i%10000==0) std::cout << i << std::endl;
		if(box.get_L()!=box2.get_L()){
			std::cout << "divergence at: " << i << std::endl;
			std::cout << "L: " << box.get_L() << std::endl;
			std::cout << "L2: " << box2.get_L() << std::endl;
			std::cout << "del: " << del << std::endl;
			std::cout << "del2: " << del2 << std::endl;
			std::cout << "del-del2: " << del-del2 << std::endl;
			for(int i=0; i<box.get_N(); i++){
				std::cout << i << ":";
				for(int j=0; j<pot2.get_cons(i).size(); j++)
					std::cout << " " << pot2.get_cons(i)[j];
				std::cout << std::endl;
			}
			break;
		}
		i++;
	}
	gsl_vector_free(r); gsl_vector_free(rsec);
	gsl_vector_free(r2); gsl_vector_free(rsec2);
	std::cout << i << " out of " << imax << std::endl;
	return i>=imax;
}

void CompareHPs::move(double alphsig, int num) {
	translate(alphsig,num);
	if(box.sym()>1){rotate(alphsig,num);}
}

void CompareHPs::translate(double alphsig, int num) {
	if(num==1){
		double L = box.get_L(); double DL = pot.get_DL();
		//double moveL;//to figure out how to move L (see 7-8-15 page 1);
		gsl_vector * mover = gsl_vector_alloc(box.get_dim());
		for(int i=0; i<box.get_N(); i++){
			gsl_vector_memcpy(mover,pot.get_DU(i));
			gsl_vector_scale(mover,alphsig);
			gsl_vector_sub(box.get_1_pos(i),mover);
			box.modify(i);//makes it all mod 1
		}
		gsl_vector_free(mover);
		//if(DL>=0) moveL=std::min(alphsig*DL,0.1*L);
		//else moveL=std::max(alphsig*DL,alphL*DL);
		box.set_L(L-alphsig*DL);
	}
	else{
		double L = box2.get_L(); double DL = pot2.get_DL();
		//double moveL;//to figure out how to move L (see 7-8-15 page 1);
		gsl_vector * mover = gsl_vector_alloc(box2.get_dim());
		for(int i=0; i<box2.get_N(); i++){
			gsl_vector_memcpy(mover,pot2.get_DU(i));
			gsl_vector_scale(mover,alphsig);
			gsl_vector_sub(box2.get_1_pos(i),mover);
			box2.modify(i);//makes it all mod 1
		}
		gsl_vector_free(mover);
		//if(DL>=0) moveL=std::min(alphsig*DL,0.1*L);
		//else moveL=std::max(alphsig*DL,alphL*DL);
		box2.set_L(L-alphsig*DL);
	}
}

void CompareHPs::rotate(double alphsig, int num) {
	if(num==1){
		gsl_vector * mover = gsl_vector_alloc(box.get_dim());
		for(int i=0; i<box.get_N(); i++){
			gsl_vector_memcpy(mover,pot.get_D_u(i));
			gsl_vector_scale(mover,alphsig);
			gsl_vector_sub(box.get_u(i),mover);
			if(box.sym()==3){
				gsl_vector_memcpy(mover,pot.get_D_v(i));
				gsl_vector_scale(mover,alphsig);
				gsl_vector_sub(box.get_v(i),mover);
			}
			box.normalize(i);
		}
		gsl_vector_free(mover);
	}
	else{
		gsl_vector * mover = gsl_vector_alloc(box2.get_dim());
		for(int i=0; i<box2.get_N(); i++){
			gsl_vector_memcpy(mover,pot2.get_D_u(i));
			gsl_vector_scale(mover,alphsig);
			gsl_vector_sub(box2.get_u(i),mover);
			if(box2.sym()==3){
				gsl_vector_memcpy(mover,pot2.get_D_v(i));
				gsl_vector_scale(mover,alphsig);
				gsl_vector_sub(box2.get_v(i),mover);
			}
			box2.normalize(i);
		}
		gsl_vector_free(mover);
	}
}

void CompareHPs::store_deriv(gsl_vector* r, int scale, int num) {
	if(num==1){
		int i=0;
		for(int j=0; j<box.get_N(); j++){
			for(int d=0; d<box.get_dim(); d++){
				gsl_vector_set(r,i,scale*gsl_vector_get(pot.get_DU(j),d));
				i++;
				if(box.sym()>1){
					gsl_vector_set(r,i,scale*gsl_vector_get(pot.get_D_u(j),d));
					i++;
					if(box.sym()==3){
						gsl_vector_set(r,i,scale*gsl_vector_get(pot.get_D_v(j),d));
						i++;
					}
				}
			}
		}
		gsl_vector_set(r,i,scale*pot.get_DL());
	}
	else{
		int i=0;
		for(int j=0; j<box2.get_N(); j++){
			for(int d=0; d<box2.get_dim(); d++){
				gsl_vector_set(r,i,scale*gsl_vector_get(pot2.get_DU(j),d));
				i++;
				if(box.sym()>1){
					gsl_vector_set(r,i,scale*gsl_vector_get(pot2.get_D_u(j),d));
					i++;
					if(box.sym()==3){
						gsl_vector_set(r,i,scale*gsl_vector_get(pot2.get_D_v(j),d));
						i++;
					}
				}
			}
		}
		gsl_vector_set(r,i,scale*pot2.get_DL());
	}
}

void CompareHPs::load_deriv(gsl_vector* r, int scale, int num) {
	if(num==1){
		int i=0;
		for(int j=0; j<box.get_N(); j++){
			for(int d=0; d<box.get_dim(); d++){
				pot.set_DU(j,d,scale*gsl_vector_get(r,i));
				i++;
				if(box.sym()>1){
					pot.set_D_u(j,d,scale*gsl_vector_get(r,i));
					i++;
					if(box.sym()==3){
						pot.set_D_v(j,d,scale*gsl_vector_get(r,i));
						i++;
					}
				}
			}
		}
		pot.set_DL(scale*gsl_vector_get(r,i));
	}
	else{
		int i=0;
		for(int j=0; j<box2.get_N(); j++){
			for(int d=0; d<box2.get_dim(); d++){
				pot2.set_DU(j,d,scale*gsl_vector_get(r,i));
				i++;
				if(box2.sym()>1){
					pot2.set_D_u(j,d,scale*gsl_vector_get(r,i));
					i++;
					if(box2.sym()==3){
						pot2.set_D_v(j,d,scale*gsl_vector_get(r,i));
						i++;
					}
				}
			}
		}
		pot2.set_DL(scale*gsl_vector_get(r,i));
	}
}

double CompareHPs::unit_fix(gsl_vector* r, int num) {
	double dotprod=0; double component;
	for(int i=0; i<r->size; i++){
		component = gsl_vector_get(r,i);
		if(i%box.sym()==0 && i<(r->size)-1) component*=box.get_L();
		dotprod+=component*component;
	}
	return dotprod;
}
