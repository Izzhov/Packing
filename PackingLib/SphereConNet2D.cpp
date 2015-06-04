/*
 * SphereConNet2D.cpp
 *
 *  Created on: Jun 2, 2015
 *      Author: izzhov
 */

#include "SphereConNet2D.h"

SphereConNet2D::SphereConNet2D(Torus2D& ptorus):torus(ptorus){
	//make con, the contact network
	for(int i=0; i<torus.get_N(); i++){
		std::vector<std::vector<int> > nex; //list of contacts for i
		con.push_back(nex); //nex will be at i
		for(int j=0; j<torus.get_N(); j++){
			if(i==j){ //can't touch itself
				if(i==(torus.get_N()-1)) break;//this means we're done
				else j++;//move on to the next one if not done
			}
			for(int k=1; k<=4; k++){
				if(torus.ell(i,j,k)<=torus.R(i,j)){//this means they're touchin
					std::vector<int> touch;
					touch.push_back(j); touch.push_back(k);
					con.at(i).push_back(touch);
				}
			}
		}
	}
	//make rs, the vector of sets of ctr-ctr vecs for ea. particle
	for(int i=0; i<torus.get_N(); i++){
		if(con.at(i).size()==0){//if particle has no contacts
			std::vector<gsl_vector*> A;
			rs.push_back(A); //this particle's vector is size zero->floater
		}
		else{//particle has contacts
			std::vector<gsl_vector*> A;//initialize vector of ctr-ctr vex
			//find&store the vec for each contact
			for(unsigned int j=0; j<con.at(i).size(); j++){
				gsl_vector * B = gsl_vector_alloc(2);//make vector
				//rel pts from i to j;hence, need negative to pt from j to i
				gsl_vector_set(B,0,-torus.rel_x(i,con[i][j][0],con[i][j][1]));
				gsl_vector_set(B,1,-torus.rel_y(i,con[i][j][0],con[i][j][1]));
				//add this to the vector of vex
				A.push_back(B);
			}
			rs.push_back(A);
		}
	}
	//remove floaters
	rm_floaters();
}

double SphereConNet2D::num_contacts() {
	int tacts = 0;
	for(int i=0; i<torus.get_N();i++){
		tacts+=con.at(i).size();
	}
	return tacts*0.5;
}

int SphereConNet2D::desired_contacts(){
	return (2*(torus.get_N()-num_floaters())-1);
}

int SphereConNet2D::num_floaters(){
	int nflt = 0;
	for(std::vector<std::vector<std::vector<int> > >::iterator it = con.begin();
			it!=con.end();++it) if((*it).size()==0) nflt++;
	return nflt;
}

std::vector<int> SphereConNet2D::which_floaters(){
	std::vector<int> whfl;
	for(int i=0;i<torus.get_N();i++) if(con[i].size()==0) whfl.push_back(i);
	return whfl;
}

void SphereConNet2D::rm_floaters(){
	while(rm_floaters_once()){}
}

bool SphereConNet2D::rm_floaters_once() {
	bool removed = false;//whether any were removed
	for(int i=0; i<torus.get_N(); i++)
		//only need to check the ones we haven't already set as floaters
		if(con.at(i).size()!=0)
			if(is_floater(i)){
				 rm(i); removed = true;
			}
	return removed;
}

bool SphereConNet2D::is_floater(int i) {
	bool fltr = false; int n = con.at(i).size();
	if(n==0||n==1) fltr = true;//always a floater when 0 or 1 contacts
	else if(n==2){//check if the two are antiparallel
		double normprod = gsl_blas_dnrm2(rs[i][0])*gsl_blas_dnrm2(rs[i][1]);
		double * dotprod; gsl_blas_ddot(rs[i][0],rs[i][1],dotprod);
		double prodsum = normprod+*dotprod;//if small->antiparallel
		//if not small->not antiparallel->is floater
		if(prodsum>std::pow(10,-14)) fltr = true;
	}
	return fltr;
}

void SphereConNet2D::rm(int i) {
	//first copy over vector of contacts to remove i from others later
	std::vector<std::vector<int> > copy = con.at(i);
	//tells us that this vector has no contacts or ctr-ctr vex any more
	std::vector<std::vector<int> > blankcon; con.at(i) = blankcon;
	std::vector<gsl_vector*> blankvec; rs.at(i) = blankvec;
	//iterate over copy to erase i from con&rs of each of its entries
	for(int j=0;j<copy.size();j++){
		int index = copy[j][0];//the index of the corr. particle
		//iterate over js contacts list to find and erase i
		for(int k=0;k<con.at(index).size();k++){
			while(con[index][k][0]==i){//must remove this
				con.at(index).erase(con.at(index).begin()+k);
				rs.at(index).erase(rs.at(index).begin()+k);
				if(k>=con.at(index).size()) break;//will segfault otherwise
			}
		}
	}
}

void SphereConNet2D::print_con(){
	for(int i=0; i<torus.get_N(); i++){
		std::cout << i << ": ";
		for(int j=0;j<con.at(i).size();j++){
			std::cout << con[i][j][0] << "," << con[i][j][1] << ";";
		}
		std::cout << "\n";
	}
}
