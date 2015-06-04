/*
 * SpheroSizes.cpp
 *
 *  Created on: May 27, 2015
 *      Author: izzhov
 */

#include "SpheroSizes.h"

SpheroSizes::SpheroSizes() {
	N = 0; totwt = 0;
}

SpheroSizes::SpheroSizes(int N, std::vector<double> rs, std::vector<double> as,
		std::vector<double> wts) {
	construct(N,rs,as,wts);
}

void SpheroSizes::construct(int N, std::vector<double> rs,
		std::vector<double> as, std::vector<double> wts) {
	set_N(N); set_rs(rs); set_as(as); set_wts(wts);
	totwt = 0;
	for(std::vector<double>::iterator it = wts.begin(); it != wts.end(); ++it)
		totwt+=*it;
	populate();
}

void SpheroSizes::populate() {
	std::vector<double> sumpts; // exact # of particles
	std::vector<int> numpts;   // # rounded down
	std::vector<double> probs;
	int totnum = 0;
	double totprob = 0;
	boost::mt19937 rng;
	unsigned int t = static_cast<unsigned int>(time(0));
	rng.seed(t);
	//calculate how many particles of each size
	for(unsigned int i=0; i < wts.size();i++){
		sumpts.push_back(N*wts.at(i)/totwt);
		numpts.push_back(floor(sumpts.at(i)));
		probs.push_back(sumpts.at(i) - numpts.at(i));
		totnum += numpts.at(i);
		totprob += probs.at(i);
	}
	//if particle fracs don't neatly divide total number of parts
	boost::random::discrete_distribution<> dist(probs);
	for(int i=0; i<N-totnum; i++){
		int q = dist(rng);
		numpts.at(q)++;
	}
	//store the radii in an array
	for(unsigned int i=0; i< wts.size(); i++)
		for(int j=0; j<numpts[i]; j++) radii.push_back(rs.at(i));
	//same with alphas
	for(unsigned int i=0; i< wts.size(); i++)
		for(int j=0; j<numpts[i]; j++) lengths.push_back(as.at(i));
}

void SpheroSizes::set_N(int N) {this->N = N;}

void SpheroSizes::set_rs(std::vector<double> rs) {this->rs = rs;}

void SpheroSizes::set_as(std::vector<double> as) {this->as = as;}

void SpheroSizes::set_wts(std::vector<double> wts) {this->wts = wts;}
