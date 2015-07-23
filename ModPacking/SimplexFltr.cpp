/*
 * SimplexFltr.cpp
 *
 *  Created on: Jul 16, 2015
 *      Author: izzhov
 */

#include "SimplexFltr.h"

SimplexFltr::SimplexFltr(std::vector<gsl_vector*> fis,
		std::vector<gsl_vector*> flocs, int d, int sym) {
	this->d=d; c=fis.size();
	int dim=fis[0]->size;
	gsl_vector * torcom = gsl_vector_alloc(3);//unit vector for below
	if(dim!=2 && sym==2){
		//then torque has 2 components->need to store first component unit vector
		torcom = mm::cross(2,flocs[0],fis[0]);
		mm::normalize(torcom);
	}

	//start by making the fti's
	if(d>dim){//then there exists torque, which you must find and paste into fti's
		for(int i=0; i<c; i++){//for each fi, make the fti
			gsl_vector * nextfti = gsl_vector_alloc(d);
			for(int j=0; j<dim; j++) //paste in the fi's
				gsl_vector_set(nextfti,j,gsl_vector_get(fis[i],j));
			if(dim==2){//then torque only has 1 component
				double newtaun=gsl_vector_get(flocs[i],0)*gsl_vector_get(fis[i],1)-
						gsl_vector_get(flocs[i],1)*gsl_vector_get(fis[i],0);//crossprod
				gsl_vector_set(nextfti,dim,newtaun);
			}
			else if(sym==2){//sym=2 in 3D->torque has 2 components
				gsl_vector * newtauv = mm::cross(2,flocs[i],fis[i]);
				//break this up into 2 components (and torque magnitude)
				double tmag; double comp1; double comp2;
				gsl_blas_ddot(newtauv,newtauv,&tmag);
				gsl_blas_ddot(newtauv,torcom,&comp1);
				if(comp1*comp1<tmag) comp2=std::sqrt(tmag-comp1*comp1);
				else comp2=0;
				gsl_vector_set(nextfti,dim,comp1);
				gsl_vector_set(nextfti,dim+1,comp2);
				gsl_vector_free(newtauv);
			}
			else{//sym=3 in 3D->torque has 3 components
				gsl_vector * newtauv = mm::cross(2,flocs[i],fis[i]);
				for(int j=dim; j<dim+newtauv->size; j++)
					gsl_vector_set(nextfti,j,gsl_vector_get(newtauv,j-dim));
				gsl_vector_free(newtauv);
			}
			ftis.push_back(nextfti);
			gsl_vector_free(fis[i]); gsl_vector_free(flocs[i]);
		}
	}
	else{//ftis are just the same as fis (no torque)
		ftis=fis;
		for(int i=0; i<c; i++) gsl_vector_free(flocs[i]);
	}
	gsl_vector_free(torcom);
	/*
	//first test: see if it makes fti's correctly
	for(int i=0; i<c; i++){
		for(int j=0; j<d; j++) std::cout << gsl_vector_get(ftis[i],j) << " ";
		std::cout << std::endl;
	}
	exit(1);
	*/
	//now initialize the gjs and bi to zeros
	for(int i=0; i<d; i++){
		bi.push_back(0);
		gsl_vector * nextgj = gsl_vector_calloc(d+c+1);
		gjs.push_back(nextgj);
	}
	//the final gj (x)
	gsl_vector * nextgj = gsl_vector_calloc(d+c+1);
	gjs.push_back(nextgj);

	//zero tolerance
	signtol = std::pow(10,-10);

	//maximum number of simplex algo steps
	maxsteps=10000;
}

bool SimplexFltr::is_floater() {
	bool isfltr=false;
	for(int i=0; i<d; i++)
		if(is_floater(i,1) || is_floater(i,-1)){
			isfltr=true; break;
		}
	//now free all of the gsl_vex that remain
	for(int i=0; i<ftis.size(); i++) gsl_vector_free(ftis[i]);
	for(int i=0; i<gjs.size(); i++) gsl_vector_free(gjs[i]);
	return isfltr;
}

bool SimplexFltr::is_floater(int unitj, int sign){
	bool isfltr=false;
	//now construct the gjs
	//j and i subscripts based on 7-16-15 notes convention:
	//j is row of tableau, i is column
	for(int j=0; j<d; j++){
		for(int i=0; i<d; i++){//first the unit left
			gsl_vector_set(gjs[j],i,j==i);
		}
		for(int i=d; i<d+c; i++){//now put in the jth entries of each fti
			double entry = gsl_vector_get(ftis[i-d],j);
			if(j==unitj) entry*=sign;//this is (P) from notes
			gsl_vector_set(gjs[j],i,entry);
		}
		//now put in the unit vector component (final column)
		gsl_vector_set(gjs[j],d+c,j==unitj);
	}
	//now make the j==d vector (aka x)
	for(int i=0; i<d; i++){//first the zeros on the left
		gsl_vector_set(gjs[d],i,0);
	}
	for(int i=d; i<d+c+1; i++){//-1*sum over j of ith entries of gjs
		double sum=0;
		for(int j=0; j<d; j++) sum+=gsl_vector_get(gjs[j],i);
		sum*=-1;
		gsl_vector_set(gjs[d],i,sum);
	}
	//make bi
	for(int j=0; j<d; j++){
		bi.at(j) = j;
	}
	/*
	//second test: see if initial tableau was made correctly
	print_tab();
	//exit(0);
	*/
	//now implement actual simplex algo part
	int indent; //the index of the entering variable
	int numsteps=0;//don't let it go over 10k
	while((indent=entering_col())<d+c){
		std::vector<std::vector<double> > rats;//ratios for simplex ratio test
		for(int j=0; j<d; j++){//store ratios
			//denominator in ratio test
			double denom=gsl_vector_get(gjs[j],indent);
			//denom must be non-zero and non-negative
			if(denom>=signtol){//store the ratio (1) and its row (0)
				std::vector<double> nextrat;
				nextrat.push_back(j);
				nextrat.push_back(gsl_vector_get(gjs[j],d+c)/denom);
				rats.push_back(nextrat);
			}
		}
		/*
		//print ratios to debug
		std::cout << "ratios: " << std::endl;
		for(int k=0; k<rats.size(); k++){
			std::cout << rats[k][0] << " " << rats[k][1] << std::endl;
		}
		*/
		if(rats.size()==0){
			std::cout << "no ratios found" << std::endl;
			exit(1);
		}
		std::vector<int> smallest;//rats-indices of smallest ratios
		smallest.push_back(0);
		//find the rat-indices of the smallest rat.s
		for(int k=1; k<rats.size(); k++){
			double diff = rats[k][1]-rats[smallest[0]][1];
			if(diff<-signtol){//neg. diff means k-th rat is smaller
				smallest.clear(); smallest.push_back(k);
			}
			else if(std::abs(diff)<signtol){//rats are equal
				smallest.push_back(k);
			}
		}
		int smlsml=0;//smlst-index with smlst basic variable row
		//find the smallest smallest (lowest basic variable col #)
		for(int k=1; k<smallest.size(); k++){
			int current = bi[rats[smallest[smlsml]][0]];
			int challenger = bi[rats[smallest[k]][0]];
			if(challenger<current) smlsml=k;
		}
		pivot(rats[smallest[smlsml]][0],indent);
		//enter entering and leave leaving
		bi.at(rats[smallest[smlsml]][0])=indent;
		numsteps++;
		if(numsteps>=maxsteps){
			std::cout << "max # of steps reached" << std::endl;
			exit(1);
		}
		//sanity check: RHS should always be positive
		for(int j=0; j<d; j++){
			double rhs=gsl_vector_get(gjs[j],d+c);
			if(rhs<=-signtol){//then neg. rhs entry exists->error
				std::cout << "negative rhs:" << std::endl;
				print_tab();
				exit(1);
			}
		}
		/*
		//third test: see if tableau is right after i iterations
		print_tab();
		//exit(1);
		*/
	}
	//bottom right corner of the tableau is the value of the obj. fxn.
	double objfxn = gsl_vector_get(gjs[d],d+c);
	if(objfxn>signtol){//objective fxn should never go positive
		std::cout << "objective function somehow positive" << std::endl;
		exit(1);
	}
	if(objfxn<-signtol){//not a floater iff objfxn is zero
		isfltr=true;
	}
	return isfltr;
}

int SimplexFltr::entering_col() {
	int entcol = 0;
	for(int i=0; i<d+c; i++){
		if(gsl_vector_get(gjs[d],i)<-signtol){
			entcol=i; break;
		}
		entcol=d+c;//will break out of while loop if can't find neg. col.
	}
	return entcol;
}

void SimplexFltr::pivot(int jj, int ii) {
	gsl_vector_scale(gjs[jj],1.0/gsl_vector_get(gjs[jj],ii));
	for(int j=0; j<d+1; j++){
		if(j==jj) j++;
		double scale = gsl_vector_get(gjs[j],ii);
		if(std::abs(scale)>=signtol){//need to subtract off this row
			gsl_vector_scale(gjs[jj],scale);
			gsl_vector_sub(gjs[j],gjs[jj]);
			gsl_vector_scale(gjs[jj],1.0/scale);
		}
	}
}

void SimplexFltr::print_tab() {
	for(int j=0; j<d+1; j++){
		for(int i=0; i<d+c+1; i++) mm::print_element(gsl_vector_get(gjs[j],i),8);
		std::cout << std::endl;
	}
	std::cout << "----------" << std::endl;
}
