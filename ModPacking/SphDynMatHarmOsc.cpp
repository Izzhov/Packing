/*
 * SphDynMatHarmOsc.cpp
 *
 *  Created on: Jun 18, 2015
 *      Author: izzhov
 */

#include "SphDynMatHarmOsc.h"

SphDynMatHarmOsc::SphDynMatHarmOsc(Torus<Sphere>& bbox,HarmPot<Torus<Sphere> >& ppot,
		SphereConNet& ccn):box(bbox),pot(ppot),cn(ccn) {
	overR.set_var(0.5);
	for(int i=0; i<2*box.get_dim(); i++){//vector has (xi,yi,...)&(xj,yj,...)
		Variable newcoord(0.0);
		coords.push_back(newcoord);
	}
	dist = f.intpow(f.minus(&coords[box.get_dim()],&coords[0]),2);
	for(int i=1; i<box.get_dim(); i++){//add all the sq's of coord dists
		dist = f.plus(f.intpow(f.minus(&coords[i+box.get_dim()],&coords[i]),2),dist);
	}
	dist = f.dpow(dist,0.5);
	dist = f.times(&overR,dist);
	dist = f.minus(f.c(1.0),dist);
	dist = f.intpow(dist,2);
	dist = f.times(f.c(0.5*pot.get_kk()),dist);

	for(int i=0; i<2*box.get_dim(); i++){//take second derivatives of all
		std::vector<Function*> nextdiff;
		for(int j=0; j<2*box.get_dim(); j++){
			nextdiff.push_back((*(*dist).diff(&coords[i])).diff(&coords[j]));
		}
		diff.push_back(nextdiff);
	}
	matsize = box.get_N()*box.get_dim(); //is a matsize x matsize matrix
	dynmat = gsl_matrix_calloc(matsize,matsize);
	calc_mat();//finds the dynamical matrix
}

void SphDynMatHarmOsc::calc_mat() {
	for(int i=0; i<box.get_N(); i++){
		for(int j=0; j<cn.num_contacts(i); j++){
			if(cn.con_elem(i,j,0)>i){//so we don't store vals twice
				//set R value
				overR.set_var(1.0/box.R(i,cn.con_elem(i,j,0)));
				//get xj,yj,(zj)
				gsl_vector * ellvec = box.ell_vec(i,cn.con_elem(i,j,0),cn.con_elem(i,j,1));
				for(int d=0; d<box.get_dim(); d++){//put these values into the last 3 coords
					coords[d+box.get_dim()].set_var(gsl_vector_get(ellvec,d));
				}
				//now we have to store the derivative values in the matrix
				for(int a=0; a<2*box.get_dim(); a++){
					for(int b=0; b<2*box.get_dim(); b++){
						int i1; int i2; //figure out which matrix indices to put in
						if(a<box.get_dim()) i1 = i*box.get_dim()+a;
						else i1 = cn.con_elem(i,j,0)*box.get_dim()+a-box.get_dim();
						if(b<box.get_dim()) i2 = i*box.get_dim()+b;
						else i2 = cn.con_elem(i,j,0)*box.get_dim()+b-box.get_dim();
						gsl_matrix_set(dynmat,i1,i2,
								gsl_matrix_get(dynmat,i1,i2)+(*diff[a][b]).eval());
						//a debugging thing
						/*
						if(a==0 && b==0 && i==0 && j==0){
							std::cout << overR.eval() << std::endl;
							for(int e=0; e<4; e++){
								std::cout << coords[e].eval() << std::endl;
							}
							std::cout << (*dist).eval() << std::endl;
							std::cout << (*diff[a][b]).eval() << std::endl;
						}
						*/
					}
				}
				gsl_vector_free(ellvec);
			}
		}
	}
}
