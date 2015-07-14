/*
 * ConNet.h
 *
 *  Created on: Jul 1, 2015
 *      Author: izzhov
 */

#ifndef CONNET_H_
#define CONNET_H_

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "Torus.h"
#include "Sphere.h"
#include "Spherocyl.h"
#include "mm.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <cstddef>
#include <vector>

//this class is made for analysis at the *end*
template <class B>//for Box
class ConNet {
private:
	B& box;
	std::vector<std::vector<std::vector<int> > > con;
	std::vector<std::vector<gsl_vector*> > rs;
public:
	ConNet(B& bbox);

	double num_contacts();//number o' contacts
	//number of contacts on element i
	int num_contacts(int i);
	//returns con[i][j][k]
	int con_elem(int i, int j, int k);

	void print_con();//prints contacts to terminal
	//makes a file with all the pairs of points for the contacts
	void output_con(std::ofstream& dstream);
};

template<class B>
inline ConNet<B>::ConNet(B& bbox):box(bbox){
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
				for(int ncon=0; ncon<box.ncon(); ncon++){
					if(std::sqrt(box.ell2(i,j,k,ncon))<=box.R(i,j)){
						//this means they're touchin
						std::vector<int> touch;
						touch.push_back(j); touch.push_back(k);
						touch.push_back(ncon);
						con.at(i).push_back(touch);
					}
				}
			}
		}
	}
	//make rs, the vector of sets of floc-floc vecs for ea. particle
	for(int i=0; i<box.get_N(); i++){
		if(con.at(i).size()==0){//if particle has no contacts
			std::vector<gsl_vector*> A;
			rs.push_back(A); //this particle's vector is size zero->floater
		}
		else{//particle has contacts
			std::vector<gsl_vector*> A;//initialize vector of floc-floc vex
			//find&store the vec for each contact
			for(unsigned int j=0; j<con.at(i).size(); j++){
				gsl_vector * G =
						box.ell_vec(i,con[i][j][0],con[i][j][1],con[i][j][2]);
				//note this vector points *away* from the particle
				A.push_back(G);//add this to the vector of vex
			}
			rs.push_back(A);
		}
	}
}

template<class B>
inline double ConNet<B>::num_contacts() {
	int tacts = 0;
	for(int i=0; i<box.get_N();i++){
		tacts+=con.at(i).size();
	}
	return tacts*0.5;
}

template<class B>
inline int ConNet<B>::num_contacts(int i) {
	return con[i].size();
}

template<class B>
inline int ConNet<B>::con_elem(int i, int j, int k) {
	return con[i][j][k];
}

template<class B>
inline void ConNet<B>::print_con() {
	for(int i=0; i<box.get_N(); i++){
		std::cout << i << ": ";
		for(int j=0;j<con.at(i).size();j++){
			std::cout << con[i][j][0] << "," << con[i][j][1] << ",";
			std::cout << con[i][j][2] << ",";
			std::cout << std::sqrt(box.ell2(i,con[i][j][0],con[i][j][1],con[i][j][2]))-
					box.R(i,con[i][j][0])<< ";";
		}
		std::cout << "\n";
	}
}

template<class B>
inline void ConNet<B>::output_con(std::ofstream& dstream) {//assumes the stream is already open to a file
	double overlap;
	double frac;//how far along the actual contact is
	int number = 1;//which number contact it is
	for(int i=0; i<box.get_N(); i++){
		for(int j=0;j<con.at(i).size();j++){
			if(con[i][j][0]>i){
				gsl_vector * firstloc = gsl_vector_alloc(box.get_dim());
				gsl_vector_memcpy(firstloc,box.get_1_pos(i));
				gsl_vector_scale(firstloc,box.get_L());
				gsl_vector_add(firstloc,box.F_loc(i,con[i][j][0],con[i][j][1],con[i][j][2]));
				gsl_vector * secondloc = box.ell_vec(i,con[i][j][0],con[i][j][1],con[i][j][2]);
				gsl_vector_add(secondloc,firstloc);
				overlap = std::sqrt(box.ell2(i,con[i][j][0],con[i][j][1],con[i][j][2]))-
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

#endif /* CONNET_H_ */
