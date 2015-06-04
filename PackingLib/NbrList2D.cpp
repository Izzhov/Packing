/*
 * NbrList2D.cpp
 *
 *  Created on: May 26, 2015
 *      Author: izzhov
 */

#include "NbrList2D.h"

NbrList2D::NbrList2D(Torus2D& ptorus):torus(ptorus){
	build_all();
}

bool NbrList2D::are_nbrs(int i, int j) {
	double d= 2*torus.R(i,j);
	return (torus.ell(i,j,1) <= d || torus.ell(i,j,2) <= d
			|| torus.ell(i,j,3) <= d || torus.ell(i,j,4) <= d);
}

bool NbrList2D::should_rebuild(){
	bool s = false; double currL = torus.get_L();
	double dl = std::abs(currL-L)*std::sqrt(2);
	double dr;
	for(int i=0; i<torus.get_N(); i++){
		dr = mymath::dist2d(locs.at(i).at(1),locs.at(i).at(2),
				torus.get_1x(i),torus.get_1y(i));
		if(dl + dr*currL >= torus.get_r(i)){
			s = true; break;
		}
	}
	return s;
}

void NbrList2D::build_nbrs() {
	nbrs.clear();
	std::vector<int> nbr;
	for(int i=0; i<torus.get_N(); i++){
		for(int j=i+1; j<torus.get_N(); j++){
			if(are_nbrs(i,j)){
				nbr.clear();
				nbr.push_back(i); nbr.push_back(j);
				nbrs.push_back(nbr);
			}
		}
	}
}

void NbrList2D::build_locs() {
	locs.clear();
	std::vector<double> loc;
	for(int i=0; i<torus.get_N(); i++){
		loc.clear();
		loc.push_back(torus.get_1x(i)); loc.push_back(torus.get_1y(i));
		locs.push_back(loc);
	}
}

void NbrList2D::build_L(){
	L = torus.get_L();
}

void NbrList2D::build_all(){
	build_nbrs(); build_locs(); build_L();
}
