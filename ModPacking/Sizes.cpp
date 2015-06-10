/*
 * Sizes.cpp
 *
 *  Created on: Jun 8, 2015
 *      Author: izzhov
 */

#include "Sizes.h"

Sizes::Sizes() {
	N=0; totwt=0;
}

Sizes::Sizes(int N, std::vector<std::vector<gsl_vector*> > szs, std::vector<double> wts) {
	set_N(N); set_szs(szs); set_wts(wts); totwt=0;
	for(std::vector<double>::iterator it = wts.begin(); it != wts.end(); ++it)
		totwt+=*it;
	populate();
}

void Sizes::populate() {
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
	//store the sizes in array fullsizes
	for(unsigned int i=0; i< wts.size(); i++)
		for(int j=0; j<numpts[i]; j++) fullsizes.push_back(szs.at(i));
}

void Sizes::set_N(int N) {
	this->N = N;
}

void Sizes::set_szs(std::vector<std::vector<gsl_vector*> > szs) {
	this->szs = szs;
}

void Sizes::set_wts(std::vector<double> wts) {
	this->wts = wts;
}
