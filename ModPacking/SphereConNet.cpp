/*
 * SphereConNet.cpp
 *
 *  Created on: Jun 11, 2015
 *      Author: izzhov
 */

#include "SphereConNet.h"

SphereConNet::SphereConNet(Torus<Sphere>& bbox, bool remove_em):box(bbox){
	//make con, the contact network
	for(int i=0; i<box.get_N(); i++){
		std::vector<std::vector<int> > nex; //list of contacts for i
		con.push_back(nex); //nex will be at i
		for(int j=0; j<box.get_N(); j++){
			if(i==j){ //can't touch itself
				if(i==(box.get_N()-1)) break;//this means we're done
				else j++;//move on to the next one if not done
			}
			for(int k=1; k<=mm::int_pow(2,box.get_dim());k++){
				if(std::sqrt(box.ell2(i,j,k,0))<=box.R(i,j)){//this means they're touchin
					std::vector<int> touch;
					touch.push_back(j); touch.push_back(k);
					con.at(i).push_back(touch);
				}
			}
		}
	}
	//make rs, the vector of sets of ctr-ctr vecs for ea. particle
	for(int i=0; i<box.get_N(); i++){
		if(con.at(i).size()==0){//if particle has no contacts
			std::vector<gsl_vector*> A;
			rs.push_back(A); //this particle's vector is size zero->floater
		}
		else{//particle has contacts
			std::vector<gsl_vector*> A;//initialize vector of ctr-ctr vex
			//find&store the vec for each contact
			for(unsigned int j=0; j<con.at(i).size(); j++){
				gsl_vector * B = box.ell_vec(i,con[i][j][0],con[i][j][1],0);
				//note this vector points *away* from the particle
				A.push_back(B);//add this to the vector of vex
			}
			rs.push_back(A);
		}
	}
	//remove floaters if you said to
	if(remove_em) rm_floaters();
}

double SphereConNet::num_contacts() {
	int tacts = 0;
	for(int i=0; i<box.get_N();i++){
		tacts+=con.at(i).size();
	}
	return tacts*0.5;
}

int SphereConNet::desired_contacts() {
	return (box.get_dim()*(box.get_N()-num_floaters())-box.get_dim()+1);
}

int SphereConNet::num_floaters() {
	int nflt = 0;
	for(std::vector<std::vector<std::vector<int> > >::iterator it = con.begin();
			it!=con.end();++it) if((*it).size()==0) nflt++;
	return nflt;
}

std::vector<int> SphereConNet::which_floaters() {
	std::vector<int> whfl;
	for(int i=0;i<box.get_N();i++) if(con[i].size()==0) whfl.push_back(i);
	return whfl;
}

void SphereConNet::rm_floaters() {
	while(rm_floaters_once()){}
}

bool SphereConNet::rm_floaters_once() {
	bool removed = false;//whether any were removed
	for(int i=0; i<box.get_N(); i++)
		//only need to check the ones we haven't already set as floaters
		if(con.at(i).size()!=0)
			if(is_floater(i)){
				rm(i); removed = true;
			}
	return removed;
}

bool SphereConNet::is_floater(int i) {
	bool fltr = false; int n = con.at(i).size();
	if(n==0||n==1) return true;//always a floater when 0 or 1 contacts
	else if(n==2){//check if the two are antiparallel
		return !mm::antiparallel(rs[i][0],rs[i][1]);//floater iff not antiparallel
	}
	else{//first check if any 2 are antiparallel
		for(int a=0; a<n-1; a++){
			for(int b=a+1; b<n; b++){//check each pair of vecs
				if(mm::antiparallel(rs[i][a],rs[i][b])) return false;
			}
		}
		//now know none are antiparallel->now it's dimension specific
		if(box.get_dim()==2){
			fltr = is_full_3_floater(i,rs[i].size()-1);
		}
		else{//dim==3
			// TODO is_full_floater 3D
		}
	}
	return fltr;
}

bool SphereConNet::is_full_3_floater(int i, int n) {
	if(n==2){
		return is_3_floater(i,0,1,2);
	}
	if(!is_full_3_floater(i,n-1)) return false;
	else for(int a=0; a<n-1; a++) for(int b=a+1; b<n; b++)
		if(!is_3_floater(i,a,b,n)) return false;
	return true;
}

bool SphereConNet::is_full_4_floater(int i, int n){
	// TODO is_full_floater 3D
	return false;
}

bool SphereConNet::is_3_floater(int x, int i, int j, int k) {
	double dotprod1 = 0; double dotprod2 = 0; //for storing dot products
	gsl_vector * crossprod1 = mm::cross(1,rs[x][i]);
	gsl_blas_ddot(crossprod1,rs[x][j],&dotprod1);
	gsl_blas_ddot(crossprod1,rs[x][k],&dotprod2);
	gsl_vector_free(crossprod1);
	if(mm::sgn(dotprod1)==mm::sgn(dotprod2)) return true;
	gsl_vector * crossprod2 = mm::cross(1,rs[x][j]);
	gsl_blas_ddot(crossprod2,rs[x][i],&dotprod1);
	gsl_blas_ddot(crossprod2,rs[x][k],&dotprod2);
	gsl_vector_free(crossprod2);
	if(mm::sgn(dotprod1)==mm::sgn(dotprod2)) return true;
	return false;
}

bool SphereConNet::is_4_floater(int x, int i, int j, int k, int l) {
	// TODO is_full_floater 3D
	return false;
}

void SphereConNet::rm(int i) {
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

int SphereConNet::num_contacts(int i){
	return con[i].size();
}

int SphereConNet::con_elem(int i, int j, int k){
	return con[i][j][k];
}

void SphereConNet::print_con(){
	for(int i=0; i<box.get_N(); i++){
		std::cout << i << ": ";
		for(int j=0;j<con.at(i).size();j++){
			std::cout << con[i][j][0] << "," << con[i][j][1] << ",";
			std::cout << std::sqrt(box.ell2(i,con[i][j][0],con[i][j][1],0))-
					box.R(i,con[i][j][0])<< ";";
		}
		std::cout << "\n";
	}
}

void SphereConNet::output_con(std::ofstream& dstream) {//assumes the stream is already open to a file
	double overlap;
	double frac;//how far along the actual contact is
	int number = 1;//which number contact it is
	for(int i=0; i<box.get_N(); i++){
		for(int j=0;j<con.at(i).size();j++){
			if(con[i][j][0]>i){
				gsl_vector * firstloc = gsl_vector_alloc(box.get_dim());
				gsl_vector_memcpy(firstloc,box.get_1_pos(i));
				gsl_vector_scale(firstloc,box.get_L());
				gsl_vector * secondloc = box.ell_vec(i,con[i][j][0],con[i][j][1],0);
				gsl_vector_add(secondloc,firstloc);
				overlap = std::sqrt(box.ell2(i,con[i][j][0],con[i][j][1],0))-
						box.R(i,con[i][j][0]);
				frac = box.get_r(i)/box.R(i,con[i][j][0]);
				for(int k=0; k<box.get_dim(); k++)
					dstream << gsl_vector_get(firstloc,k) << " ";
				for(int k=0; k<box.get_dim(); k++)
					dstream << gsl_vector_get(secondloc,k) << " ";
				dstream << overlap << " " << frac << " ";
				dstream << number << "\n"; number++;
				gsl_vector_free(firstloc); gsl_vector_free(secondloc);
			}
		}
	}
}
